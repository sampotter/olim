function [u, lam, deg] = F0_exact(u0, du, p0, dP, sh)
    [Q, R] = qr(dP, 0);
    I = eye(size(Q, 1));
    q = inv(R')*du/sh;
    if dot(q, q) >= 1
        deg = true;
        return;
    end
    lopt = sqrt(p0'*(I - Q*Q')*p0/(1 - dot(q, q)));
    lam = -inv(R)*(Q'*p0 + lopt*q);
    u = u0 + sh*p0'*(p0 + dP*lam)/lopt;
end
