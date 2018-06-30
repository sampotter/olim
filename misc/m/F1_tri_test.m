clear;

h = rand;
u0 = rand;
u1 = rand;
s = rand;
s0 = rand;
s1 = rand;
theta = rand;
p0 = randn(3, 1);
p1 = randn(3, 1);

fprintf(['u0 = %g, u1 = %g, h = %g, s = %g, s0 = %g, s1 = %g, theta ' ...
         '= %g\n'], u0, u1, h, s, s0, s1, theta);

fprintf('p0 = (%g, %g, %g)\n', p0(1), p0(2), p0(3));
fprintf('p1 = (%g, %g, %g)\n', p1(1), p1(2), p1(3));

[F1opt, lamopt, lam, F1iters] = F1_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1);

fprintf('F1 = %0.16g\n', F1opt);

u = @(lam) (1 - lam)*u0 + lam*u1;
s = @(lam) (1 - theta)*s + theta*((1 - lam)*s0 + lam*s1);
p = @(lam) p0 + lam*(p1 - p0);
F1 = @(lam) u(lam) + h*s(lam)*sqrt(p(lam)'*p(lam));

Lams = linspace(-0.5, 1.5, 151);
F1s = arrayfun(F1, Lams);

figure;
hold on;
plot(Lams, F1s);
scatter(lam, F1iters, 'r');
plot(lamopt, F1(lamopt), '*');
xlim([min(Lams) max(Lams)]);
