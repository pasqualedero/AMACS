function k0 = compk0(parameters)
T=parameters.T;
for k=1:100
    b=1;
    for j=1:size(T,1)
        if b==1
            b=b && jk(j,k,parameters);
        end
    end
    if b==1
        break
    end
end
if b==1
    k0=k;
else
    k0=-k;
end
end


function b=jk(j,k,parameters)

b=0;

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

%%
% w \in \W^d
con=[];
for jj=1:size(T,1)
    con=T(jj,:)*A*w<=g(jj)-d;
end

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

Rx=Rx+Phi^(k)*G;
% Solve the optimization problem
solution = optimize(con, -T(j,:)*(Hc*(Phi^k*x+Rx*w)+L*w), options);

if solution.problem~=0
    disp(solution);
    disp([j,k])
    error('error')
end

Gkj=T(j,:)*(Hc*(Phi^k*value(x)+Rx*value(w))+L*value(w))-g(j);

if Gkj<=0
    b=1;
end

end

