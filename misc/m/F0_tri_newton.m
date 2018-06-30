function [F0opt, lamopt, lam1, lam2, degenerate] = F0_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1, silent)
    if nargin < 10
        silent = true;
    end
    
    du = u1 - u0;
    dp = p1 - p0;
    savg = @(lam) (1 - lam)*s0 + lam*s1;
    stheta = @(lam) (1 - theta)*s + theta*savg(lam);
    alpha = -du/(stheta(0.5)*h);

    u = @(lam) (1 - lam)*u0 + lam*u1;
    p = @(lam) p0 + lam*dp;
    q = @(lam) p(lam)'*p(lam);
    l = @(lam) sqrt(q(lam));

    F0 = @(lam) u(lam) + h*stheta(lam)*l(lam);

    tmp = alpha^2 - norm(dp, 2)^2;
    a = norm(dp, 2)^2*tmp;
    b = dp'*p0*tmp;
    c = alpha^2*norm(p0, 2)^2 - (dp'*p0)^2;

    degenerate = false;
    if b*b < a*c
        degenerate = true;
        if F0(0) <= F0(1)
            lamopt = 0;
        else
            lamopt = 1;
        end
        if ~silent
            fprintf('degenerate problem\n');
            fprintf('lamopt = %g\n', lamopt);
        end
    else
        lhs = -b/a;
        rhs = sqrt(b^2 - a*c)/a;

        lam1 = lhs - rhs;
        lam2 = lhs + rhs;

        check = @(lam) abs(alpha*l(lam) - dp'*p(lam));
        tol = 1e-10;
        if check(lam1) < check(lam2)
            lamopt = lam1;
        else
            lamopt = lam2;
        end
        if lamopt < 0 || 1 < lamopt
            if F0(0) <= F0(1)
                lamopt = 0;
            else
                lamopt = 1;
            end
        end
        if ~silent
            fprintf('check1 = %g, check2 = %g\n', check(lam1), check(lam2));
            fprintf('lam1 = %g, lam2 = %g, lamopt = %g\n', lam1, lam2, lamopt);
        end
    end

    if ~exist('lam1')
        lam1 = NaN;
    end
    if ~exist('lam2')
        lam2 = NaN;
    end

    F0opt = F0(lamopt);
end
