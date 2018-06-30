function xopt = qpe(G, c, A, b)
% Solve an equality-constrained quadratic program using the Schur
% complement method (following Nocedal & Wright)
    tmp1 = G\(A');
    tmp2 = G\c;
    mu = (A*tmp1)\(A*tmp2);
    xopt = tmp1*mu - tmp2;
end
