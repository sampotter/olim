function [xopt, p] = qpe(G, c, A, b)
% Solve an equality-constrained quadratic program using the Schur
% complement method (following Nocedal & Wright)
    mu = pinv(A*inv(G)*A')*A*inv(G)*c;
    xopt = inv(G)*(A'*mu - c);
end
