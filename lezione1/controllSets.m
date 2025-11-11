%% Compute one-step ahead Controllable Sets
clear; close all; clc;

% Define the polytopic description

yalmip('clear'); 

parameters.Ts=0.1;
parameters.Kmin = 1.8;
parameters.Kmax = 2.8;
parameters.M=1;
parameters.b=1;

Ts=parameters.Ts;
Kmin=parameters.Kmin;
Kmax = parameters.Kmax;
M = parameters.M;
b=parameters.b;

A1=[1 Ts; -Kmin*Ts/M 1-b*Ts/M];
A2=[1 Ts; -Kmax*Ts/M 1-b*Ts/M];
A = {A1, A2};

B1=[0; Ts/M];
B = {B1, B1};

C = [1 0; 0 1];

[parameters.n,parameters.m]=size(B{1});

% Compute Epsilon0's shaping matrix, keeping in mind we impose euclidean
% norm input constraints (specified in function's arguments)
x0 = [0.1; -0.2];
[F, Qe] = constrainedRHC(A,B,C,x0,1,1);

% Create ellipsoid
figure(); hold on;

Q = inv(Qe);
el0 = Ellipsoid(Qe);
el0.plotEll(Color='r');

% Call one-step-ahead function and store shaping matrices
els = cell(1,20);
els{1} = el0;

P = cell(1,20);
P{1} = Qe;

for i=2:20
    [els{i}, flag] = Osa(els{i-1},A,B,C,1,1);
    els{i}.plotEll('Color','b');
    P{i} = els{i}.Q;
end

% Online
terminalSetInfo.inTerminalSet = 0;
terminalSetInfo.centeredInZero = 0;

x = [-0.1; 0.95];
X = zeros(parameters.n, length(P));
X(:,1) = x;
indice=2;

U = zeros(parameters.m, length(P));

u0 = 0; % Initial guess for the input u (scalar)
options = optimoptions('fmincon', 'Display', 'none'); % Hide fmincon output

disp('--- Starting Online Simulation ---');
fprintf('Initial state: x = [%.3f; %.3f]\n', x(1), x(2));

while ~ terminalSetInfo.inTerminalSet
    
    index = 0; 
    
    % Find the smallest set 'i' that contains the current state 'x'
    for i=1:length(P)
        el = Ellipsoid(P{i});
        if el.isInternal(x)
            index = i;
            break; % Found the smallest set, stop searching
        end
    end

    % Case 1: State is in the terminal set (P{1}) -> break and apply
    % feedback
    if index == 1
        terminalSetInfo.inTerminalSet = 1;
        terminalSetInfo.firstStateReached = x;
        terminalSetInfo.firstStateReachedIndex = indice - 1;
        disp('State has reached the terminal set P{1}.');
        break; 
    end
    
    % Case 2: State is outside all controllable sets
    if index == 0
        disp('Error: State is outside all controllable sets');
        break; % Stop the while loop
    end
    
    % Case 3: State is in set P{i} (where i > 1) -> compute control
    fprintf('State is in set P{%d}. Computing control to reach P{%d}.\n', index, index-1);
    
    P_target = P{index-1}; % Target the next smaller set
    
    % Create function handles with current state 'x'
    objFun = @(u_in) objective(u_in, A, B, x, P_target);
    conFun = @(u_in) mycon(u_in, 1, 1, A, B, C, x, P_target); % Using u_max=1, y_max=1
    
    [u_optimal, fval] = fmincon(objFun, u0, ...
                          [], [], ... 
                          [], [], ... 
                          [], [], ... 
                          conFun, ...  
                          options);
    
    U(:,indice) = u_optimal;
    % Apply u(t) to the system
    x = A{1}*x + B{1}*u_optimal;
    % Add new computed state to X
    X(:,indice) = x; 

    fprintf('State updated from [%.3f; %.3f] to [%.3f; %.3f] (u=%.4f)\n', ...
            X(1,indice-1), X(2,indice-1), X(1,indice), X(2,indice), u_optimal); 

    indice = indice+1;
end

% Terminal set -> feedback

while ~ terminalSetInfo.centeredInZero
    u = F*X(:,indice-1);
    x = A{1}*X(:,indice-1) + B{1}*u;
    X(:,indice) = x; 
    if vecnorm(x) <= 1e-4
        terminalSetInfo.centeredInZero = 1;
    else
        indice = indice + 1;
    end
end


plot(X(1,:),X(2,:),'Marker','+');
hold off;
