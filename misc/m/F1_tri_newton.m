function [F1opt, lamopt, lam, F1iters] = F1_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1)
    du = u1 - u0;
    dp = p1 - p0;
    ds = s1 - s0;
    s = @(lam) (1 - theta)*s + theta*((1 - lam)*s0 + lam*s1);

    p = @(lam) p0 + lam*dp;

    u = @(lam) (1 - lam)*u0 + lam*u1;

    q = @(lam) p(lam)'*p(lam);
    dq = @(lam) 2*dp'*p(lam);
    d2q = 2*dp'*dp;

    l = @(lam) sqrt(q(lam));

    F1 = @(lam) u(lam) + h*s(lam)*sqrt(p(lam)'*p(lam));
    dF1 = @(lam) du + h*(ds*theta*p(lam)'*p(lam) + s(lam)*dp'*p(lam))/l(lam);
    d2F1 = @(lam) h*(s(lam)*(dp'*p(lam))^2/(p(lam)'*p(lam)) + ds*theta*dp'*p(lam) + ...
                     2*s(lam)*dp'*dp)/(l(lam));

    g = @(lam) -dF1(lam)/d2F1(lam);

    iter = 1;
    maxiters = 100;

    lam(iter) = 0.5;
    F1iters(iter) = F1(lam(iter));
    alpha(iter) = 1;

    c = 0.4;
    c1 = 1e-4;
    c2 = 0.9;

    suff_dec_cond = @(x, a) F1(x + a*g(x)) <= F1(x) + c1*a*dF1(x)*g(x);

    goldstein_conds = @(x, a) ...
        F1(lam) + (1 - c)*a*dF1(x)*g(x) <= F1(x + a*g(x)) && ...
        F1(x + a*g(x)) <= F1(x) + c*a*dF1(x)*g(x);

    wolfe_conds = @(x, a) ...
        F1(x + a*g(x)) <= F1(x) + c1*a*dF1(x)*g(x) && ...
        dF1(x + a*g(x)) >= c2*dF1(x)*g(x);

    strong_wolfe_conds = @(x, a) ...
        F1(x + a*g(x)) <= F1(x) + c1*a*dF1(x)*g(x) && ...
        abs(dF1(x + a*g(x))) >= c2*abs(dF1(x)*g(x));

    while iter < 100 && norm(g(lam(iter)), 'inf') > eps
        alpha(iter) = 1;
        while ~suff_dec_cond(lam(iter), alpha(iter))
            alpha(iter) = 0.9*alpha(iter);
        end

        gs(iter) = g(lam(iter));
        lam(iter + 1) = lam(iter) + alpha(iter)*gs(iter);
        iter = iter + 1;
        lam(iter) = max(0, min(1, lam(iter)));

        F1iters(iter) = F1(lam(iter));
        assert(F1iters(iter) <= F1iters(iter - 1));

        tol = 10*eps;
        dlam = lam(iter) - lam(iter - 1);
        dF1iters = F1iters(iter) - F1iters(iter - 1);
        if abs(dlam) < tol || abs(dF1iters) < tol
            break
        end
    end
    lamopt = max(0, min(1, lam(iter)));
    F1opt = F1(lamopt);
end
