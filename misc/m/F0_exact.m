function [u, lam, deg] = F0_exact(u0, du, p0, dP, sh)
    deg = false;
    [Q, R] = qr(dP, 0);
    I = eye(size(Q, 1));
    q = inv(R')*du/sh;
    if dot(q, q) >= 1
        deg = true;
        u = inf;
        lam = [nan nan];
    else
        lopt = sqrt(p0'*(I - Q*Q')*p0/(1 - dot(q, q)));
        lam = -inv(R)*(Q'*p0 + lopt*q);
        if ~in_simplex(lam)
            deg = true;
            u = inf;
        else
            u0 + du'*lam + sh*lopt
            u = u0 + sh*p0'*(p0 + dP*lam)/lopt;
        end
    end
end

function [tf] = in_simplex(lam)
    tf = all(lam >= 0) & sum(lam) <= 1;
end
