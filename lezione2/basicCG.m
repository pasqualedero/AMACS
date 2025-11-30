clear; clc; close all;

%% Define system

parameters.K=1;
parameters.M=1;
parameters.b=1;
parameters.Ts=0.5;
K=parameters.K;
M=parameters.M;
b=parameters.b;
Ts=parameters.Ts;

A= [1 Ts; -K*Ts/M 1-b*Ts/M];
B= [0; Ts/M];
C = [1 0];
D = 0;
[parameters.n,parameters.m]=size(B);
parameters.S=eye(parameters.n);

%% 1 - Asymptotical Stab.
spectr = eig(A);

%% 2 - Offset free
offset = C * inv(eye(parameters.n)-A) * B;
if offset == eye(parameters.m)
    disp('satisfied');
end

%% Step response (without CG)

% Create custom signal: Simul. Time N = 200s
N = 200;
phase_duration = N/4;
steps_per_phase = phase_duration/Ts;

s1 = 0.1 * ones(1, steps_per_phase);
s2 = 0.5 * ones(1, steps_per_phase);
s3 = -0.5 * ones(1, steps_per_phase);
s4 = 0.8 * ones(1, steps_per_phase);

signal = [s1, s2, s3, s4];

% Total steps
total_steps = length(signal);

% Step indexes: from 0 to N with Ts as increment
k = 0 : Ts : (total_steps - 1) * Ts;

% Simulate
sys = ss(A,B,C,D,Ts);
[y,~,x] = lsim(sys,signal',k);

figure('Name','Input Signal vs Input Response')
hold on
plot(k,signal,'LineStyle','-.','LineWidth',1)
plot(k,y,'LineStyle','-','LineWidth',1)
grid on
legend('Input Signal','Input Response')
hold off

%% ---------- CG Application (no disturbances) --------------------------------------------

%% OFFLINE PHASE
% Define matrices

% constraints matrix T
parameters.T = [1 0 0;
                -1 0 0;
                0 1 0;
                0 -1 0;
                0 0 1;
                0 0 -1];
parameters.g = ones(6,1);
parameters.Hc = [1 0; 0 1; 0 0];    % related to x(t)
parameters.Phi = A;
parameters.G = B;
parameters.L = [0; 0; 1];           % related to g(t)
parameters.delta = 0.1;

% Compute k0

k0 = compk0(parameters);

%% ONLINE PHASE

% matrix \Psi
parameters.Psi = 1;

% store x and y and g
X = zeros(parameters.n, length(signal))';
Y = zeros(parameters.m, length(signal))';
g_store = zeros(parameters.m, length(signal));

X(1,:) = [0 0];

% run simulation
for i = 1:length(signal)-1
    g = computeG(signal(i), X(i,:)', k0, parameters);

    x_next = A * X(i,:)' + B * g';
    y_next = C * X(i,:)' + D * g';

    X(i+1,:) = x_next';
    Y(i+1,:) = y_next';

    g_store(:,i) = g;

end

%% PLOTS

% Non-CG Approach

figure('Name','CG vs Non-CG controlled system');
hold on

subplot(6,3,[1 2 4 5 7 8])
hold on;
plot(k,signal,'LineStyle','--','LineWidth',1.5,'Color','r');
plot(k,y,'LineWidth',2,'Color','c');
grid on;
title('Non-CG approach')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
legend('reference','response');
ylim([-2 2]);
hold off

subplot(6,3,3)
plot(k,x(:,1),'LineWidth',1.5,'Color','b')
grid on;
title('x_1')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

subplot(6,3,6)
plot(k,x(:,2),'LineWidth',1.5,'Color','b')
grid on;
title('x_2')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

subplot(6,3,9);
plot(k,signal,'LineWidth',1.5,'Color','b')
grid on;
title('u(t)')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

% CG Approach

subplot(6,3,[10 11 13 14 16 17])
hold on;
plot(k,signal,'LineStyle','--','LineWidth',1.5,'Color','r');
plot(k,Y,'LineWidth',2,'Color','c');
grid on;
title('CG approach')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
legend('reference','response');
ylim([-2 2]);
hold off

subplot(6,3,12)
plot(k,X(:,1),'LineWidth',1.5,'Color','b')
grid on;
title('x_1')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

subplot(6,3,15)
plot(k,X(:,2),'LineWidth',1.5,'Color','b')
grid on;
title('x_2')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

subplot(6,3,18);
plot(k,g_store,'LineWidth',1.5,'Color','b')
grid on;
title('u(t)')
xlabel('k-step')
yline(1,'LineStyle','--','Color','g','LineWidth',1.5);
yline(-1,'LineStyle','--','Color','g','LineWidth',1.5);
ylim([-2 2]);

hold off


