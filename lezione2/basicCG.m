clear; clc;

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

% Create custom signal
N = 1000;
signal = zeros(1,N);
k = 1:Ts:N;

signal(k<200) = 10;
signal(k>=200 & k<500) = 20;
signal(k>=500 & k<=700) = 5;
signal(k>700) = 15;

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

% run simulation
for i = 1:Ts:N
    g = computeG(signal(i), x(i,:)', k0, parameters);
end


