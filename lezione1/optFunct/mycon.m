function [c,ceq] = mycon(x,u_max,y_max,A,B,C,xt,P1)
    n = length(A);
    c = zeros(n+2,1);

    xt_evol_max = zeros(size(xt));
    xt_evol_max_norm = -inf;
    
    % dynamic constraints
    for i=1:n
        row = (A{i}*xt + B{i}*x)'*P1*(A{i}*xt + B{i}*x)-1;
        c(i,1) = row;

        xt_evol = (A{i}*xt + B{i}*x);
        xt_evol_norm = vecnorm(xt_evol);

        if xt_evol_norm > xt_evol_max_norm
            xt_evol_max_norm = xt_evol_norm;
            xt_evol_max = xt_evol;
        end
        
    end
    % input constraint
    c(n+1,1)= vecnorm(x)-u_max;
    % output constraint 
    c(n+2) = vecnorm(C * xt_evol_max) - y_max;
    ceq = [];
end