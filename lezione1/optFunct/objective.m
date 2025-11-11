function [f] = objective(x,A,B,xt,Q1)
% max (Aj * xt + Bj * x)' * Q1 * (Aj * xt + Bj * x)
    n = length(A);
    f_max = zeros(n,1);
    for i=1:n
        row = (A{i}*xt+B{i}*x)'*Q1*(A{i}*xt+B{i}*x);
        f_max(i,1) = row;
    end
    f = max(f_max);
end