clear;

quad = 'mp0';
fac = true;
bv = false;
test_sqp = true;
tol = eps;

if fac, assert(~bv); end

I = eye(3);

h = 1;

u = sqrt(2);
u0 = 0;
u1 = 1;
u2 = sqrt(2);
du = [u1 - u0; u2 - u0];

s = 1;
s0 = 1;
s1 = 1;
s2 = 1;
ds = [s1 - s0; s2 - s0];
sf = 1;

% P = rand(3);
P = [1 1 1; 0 0 1; -1 0 0];
% P = I;
% P = [1 1 1; 0 1 1; 0 0 1];
% P = [1 1 0; 1 0 1; 0 1 1];
if bv, assert(all(P(:) == 0 | P(:) == 1)); end
assert(rank(P) == 3);

% pf = randn(3, 1);
pf = [1; 0; -1];

p0 = P(:, 1);
dP = P(:, 2:3) - p0;
p = @(x) p0 + dP*x;
l = @(x) norm(p(x));
n = @(x) p(x)/l(x);
lf = @(x) norm(p(x) - pf);
nf = @(x) ((p(x) - pf)/max(eps, lf(x)))*(lf(x) > 1e1*eps);
C = @(x) I - n(x)*n(x)';
Cf = @(x) I - nf(x)*nf(x)';

T = @(x) sf*h*lf(x);
tau0 = u0 - T([0; 0]);
tau1 = u1 - T([1; 0]);
tau2 = u2 - T([0; 1]);
dtau = [tau1 - tau0; tau2 - tau0];
tau = @(x) tau0 + dtau'*x;

if strcmp(quad, 'rhr')
    if fac
        f = @(x) tau(x) + T(x) + s*h*l(x);
        df = @(x) dtau + sf*h*dP'*nf(x) + s*h*dP'*n(x);
        d2f = @(x) s*h*dP'*C(x)*dP/l(x) + ...
              (sf*h*dP'*Cf(x)*dP/max(eps, lf(x)))*(lf(x) > 1e1*eps);
    else
        f = @(x) u0 + du'*x + s*h*l(x);
        df = @(x) du + s*h*dP'*n(x);
        d2f = @(x) s*h*dP'*C(x)*dP/l(x);
    end
elseif strcmp(quad, 'mp0')
    sh = h*(s + (s0 + s1 + s2)/3)/2;
    if fac
        f = @(x) tau(x) + T(x) + sh*l(x);
        df = @(x) dtau + sf*h*dP'*nf(x) + sh*dP'*n(x);
        d2f = @(x) sh*dP'*C(x)*dP/l(x) + ...
              (sf*h*dP'*Cf(x)*dP/max(eps, lf(x)))*(lf(x) > 1e1*eps);
    else
        f = @(x) u0 + du'*x + sh*l(x);
        df = @(x) du + sh*dP'*n(x);
        d2f = @(x) sh*dP'*C(x)*dP/l(x);
    end
elseif strcmp(quad, 'mp1')
    sh = @(x) (s + s0 + ds'*x)*h/2;
    if fac
        f = @(x) tau(x) + T(x) + sh(x)*l(x);
        df = @(x) dtau + sf*h*dP'*nf(x) + l(x)*ds*h/2 + sh(x)*dP'*n(x);
        d2f = @(x) (sf*h*dP'*Cf(x)*dP/max(eps, lf(x)))*(lf(x) > 1e1*eps) + ...
              (h*(dP'*p(x)*ds' + ds*p(x)'*dP)/2 + sh(x)*dP'*C(x)*dP)/l(x);
    else
        f = @(x) u0 + du'*x + sh(x)*l(x);
        df = @(x) du + l(x)*ds*h/2 + sh(x)*dP'*n(x);
        d2f = @(x) (h*(dP'*p(x)*ds' + ds*p(x)'*dP)/2 + sh(x)*dP'*C(x)*dP)/l(x);
    end
else
    error('not implemented yet');
end

fprintf('u0 = %0.16g;\n', u0);
fprintf('u1 = %0.16g;\n', u1);
fprintf('u2 = %0.16g;\n', u2);
fprintf('h = %0.16g;\n', h);
fprintf('s = %0.16g;\n', s);
fprintf('s0 = %0.16g;\n', s0);
fprintf('s1 = %0.16g;\n', s1);
fprintf('s2 = %0.16g;\n', s2);
if fac
    fprintf('sf = %0.16g;\n', sf);
end
if ~bv
    for i = 1:3
        for j = 1:3
            fprintf('p%d[%d] = %0.16g;\n', i - 1, j - 1, P(j, i));
        end
    end
end
if fac
    for i = 1:3
        fprintf('pf[%d] = %0.16g;\n', i - 1, pf(i));
    end
end
if ~test_sqp
    fprintf('lam[0] = %0.16g;\n', lam(1));
    fprintf('lam[1] = %0.16g;\n', lam(2));
end
if fac
    fprintf('F_fac_wkspc<%s, 2> w;\n', upper(quad));
else
    fprintf('F_wkspc<%s, 2> w;\n', upper(quad));
end
if fac
    fprintf(['set_args<%s, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, ' ...
             's2, h, pf, sf);\n'], upper(quad));
else
    if bv
        fprintf(['set_args<%s, 3, P%d%d%d, P%d%d%d, P%d%d%d>(w, u0, u1, u2, ' ...
                 's, s0, s1, s2, h);\n'], upper(quad), P(:, 1), P(:, 2), ...
                P(:, 3));
    else
        fprintf(['set_args<%s, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, ' ...
                 's2, h);\n'], upper(quad));
    end
end
if test_sqp
    lam0 = [1; 1]/3;
    A = [eye(2); -ones(1, 2)];
    b = [zeros(2, 1); -1];
    out = sqp(f, df, d2f, A, b, lam0, eps, 100);
    if bv
        fprintf(['cost_functor_bv<%s, 3, P%d%d%d, P%d%d%d, P%d%d%d> ' ...
                 'func {w};\n'], upper(quad), P(:, 1), P(:, 2), P(:, 3));
    else
        if fac
            fprintf('cost_functor_fac<%s, 3> func {w, p0, p1, p2, pf};\n', ...
                    upper(quad));
        else
            fprintf('cost_functor<%s, 3> func {w, p0, p1, p2};\n', ...
                    upper(quad));
        end
    end
    fprintf('updates::info<2> info;\n');
    fprintf('sqp_bary<decltype(func), 3, 2>()(func, info.lambda, &error);\n');
    fprintf('ASSERT_FALSE(error);\n');
    fprintf('ASSERT_NEAR(info.lambda[0], %0.16g, %g);\n', out.xs(1, out.iters + 1), tol);
    fprintf('ASSERT_NEAR(info.lambda[1], %0.16g, %g);\n', out.xs(2, out.iters + 1), tol);
    if exist('u')
        fprintf('ASSERT_NEAR(info.value, %0.16g, %g);\n', out.fs(out.iters  + 1), tol);
    end
else
    lam = rand(2, 1);
    lam = lam/sum(lam);
    f_lam = f(lam);
    df_lam = df(lam);
    d2f_lam = d2f(lam);
    if fac
        fprintf('set_lambda<%s, 3>(w, p0, p1, p2, pf, lam);\n', ...
                upper(quad));
    else
        if bv
            fprintf('set_lambda<%s, 3, P%d%d%d, P%d%d%d, P%d%d%d>(w, lam);\n', ...
                    upper(quad), P(:, 1), P(:, 2), P(:, 3));
        else
            fprintf('set_lambda<%s, 3>(w, p0, p1, p2, lam);\n', ...
                    upper(quad));
        end
    end
    fprintf('eval(w, f);\n');
    fprintf('grad(w, df);\n');
    fprintf('hess(w, d2f);\n');
    fprintf('ASSERT_NEAR(f, %0.16g, %g);\n', f_lam, tol);
    fprintf('ASSERT_NEAR(df[0], %0.16g, %g);\n', df_lam(1), tol);
    fprintf('ASSERT_NEAR(df[1], %0.16g, %g);\n', df_lam(2), tol);
    fprintf('ASSERT_NEAR(d2f[0], %0.16g, %g);\n', d2f_lam(1), tol);
    fprintf('ASSERT_NEAR(d2f[1], %0.16g, %g);\n', d2f_lam(2), tol);
    fprintf('ASSERT_NEAR(d2f[2], %0.16g, %g);\n', d2f_lam(4), tol);
end