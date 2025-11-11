function [ep,s] = Osa(e0,A,B,C,umax,ymax)
%OSA Summary of this function goes here
%   Given an ell e0, this function computes the inner ellipsoidal 
% approximation of its one-step-ahed controllable set
[n,m]=size(B{1});
P0=e0.Q;
Qvar=sdpvar(n+m);
con=Qvar>=0;

if isempty(C)
    C=eye(n);
end

for j=1:length(A)
    P1ext_j=[A{j}';B{j}']*P0*[A{j} B{j}];
    [U,Sj,~]=svd(P1ext_j);
    Q1ext_j=pinv(P1ext_j);
    r=rank(P1ext_j);
    Lat=[eye(r) zeros(r,n+m-r)];
    con=[con Lat*U'*Qvar*U*Lat'<=inv(Sj(1:r,1:r))];
end

P1ext_vin=[C'*C/(ymax^2) zeros(n,m); zeros(m,n) eye(m)/(umax^2)];
%P1ext_vin=[zeros(n) zeros(n,m); zeros(m,n) eye(m)/(umax^2)];
[U,S,~]=svd(P1ext_vin);
Q1ext_vin=pinv(P1ext_vin);
r=rank(P1ext_vin);
Lat=[eye(r) zeros(r,n+m-r)];
con=[con Lat*U'*Qvar*U*Lat'<=inv(S(1:r,1:r))];

% Nesting
% P1ext_nes=[P0 zeros(n,m); zeros(m,n) zeros(m)];
% [U,S,~]=svd(P1ext_nes);
% Q1ext_nes=pinv(P1ext_nes);
% r=rank(P1ext_nes);
% Lat=[eye(r) eye(r,n+m-r)];
% con=[con Lat*U*Qvar*U'*Lat'>=inv(S(1:r,1:r))];
con=[con Qvar(1:n,1:n)>=pinv(P0)];

% Define solver options 
options = sdpsettings( 'verbose', 0);

% Solve the optimization problem
sol = optimize(con, -log(det(Qvar(1:n,1:n))), options);

% Check if the problem is solved
if sol.problem == 0
    Qe = value(Qvar);    
    s=0;
    ep=Ellipsoid(inv(Qe(1:n,1:n)));
else
    disp(sol.info);
    ep=e0;
    s=1;
    ep=[];
    disp('Something went wrong:');
end



end

