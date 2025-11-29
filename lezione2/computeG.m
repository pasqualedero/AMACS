function g = computeG(reference,x,k0,parameters)
% returns the command g to apply to primal controller

T=parameters.T;
g=parameters.g;
Hc=parameters.Hc;
Phi=parameters.Phi;
G=parameters.G;
L=parameters.L;
d=parameters.delta;
Psi = parameters.Psi;

Rx = zeros(parameters.n,parameters.m);

yalmip('clear');
w = sdpvar(parameters.m,1);

con = T * (Hc * inv(eye(size(Phi)) - Phi) * G + L <= g - d * ones(size(g)));
for k = 0:k0
    Rc=Hc*Rx+L;
    con = [con T*Hc*Phi^k*x+T*Rc*w<=g];
    Rx=Rx+Phi^k*G;
end

options = sdpsettings('solver','sedumi','verbose', 0);
solution = optimize(con, (w-reference)'*Psi(w-reference), options);

if solution.problem~=0
    disp(solution);
    error('error')
end

end