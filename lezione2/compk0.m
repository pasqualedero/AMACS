function k0 = compk0(parameters)
T=parameters.T;
k0 = -1;
for k=1:100
    moveOn = 1;
    for j = 1:size(T,1)
        if moveOn==1
            moveOn = Gjk(j,k,parameters);
        else
            break;
        end
    end
    if moveOn
        k0 = k;
        break;
    end
end
end


function lt0 = Gjk(j,k,parameters)
% Compute the optimizazion problem, where j,k is fixed.
% Impose:
%   - 1 constraint for set membership of w
%   - k-1 contraints related to the fact that constraints were satisfied
%     in the past from 0 to k-1
%
% INPUTS
% j     j-th constraint    
% k     k step
% OUTPUTS
% lt0   Gjk <= 0 : boolean

lt0=0;

T=parameters.T;
g=parameters.g;
Hc=parameters.Hc;
Phi=parameters.Phi;
G=parameters.G;
L=parameters.L;
d=parameters.delta;

A=Hc*inv(eye(parameters.n)-Phi)*G+L;

yalmip('clear')
x=sdpvar(parameters.n,1);
w=sdpvar(parameters.m,1);

% CONSTRAINTS

% subject to: w in W^d
con = T*A*w <= g-d*ones(size(g));
% con = []; 
% for jj=1:size(T,1)
%     con= T(jj,:)*A*w<=g(jj)-d;
% end

% subject to: all constraints were prev. satisfied up to k-1 included 
for jj = 1:size(T,1)
    Rx=zeros(parameters.n,parameters.m);
    for i=0:k-1
        Rc=Hc*Rx+L;
        %xk=(Phi^k*x+Rx*w);
        %ck=Hc*xk+L*w;
        con=[con T(jj,:)*(Hc*Phi^i*x+Rc*w)<=g(jj)];
        Rx=Rx+Phi^i*G;
    end
end

% Define solver options
options = sdpsettings('solver','sedumi','verbose', 0);

%con=[con abs(x(1))<=parameters.xmax abs(x(2))<=parameters.vmax abs(w)<=parameters.umax]

%Rx=Rx+Phi^(k)*G;
% Solve the optimization problem
solution = optimize(con, -T(j,:)*(Hc*(Phi^k*x+Rx*w)+L*w)-g(j), options);

if solution.problem~=0
    disp(solution);
    disp([j,k])
    error('error')
end

Gkj=T(j,:)*(Hc*(Phi^k*value(x)+Rx*value(w))+L*value(w))-g(j);

if Gkj<=0
    lt0 = 1;
end

end


