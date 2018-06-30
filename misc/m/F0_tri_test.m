clear;

%%% TODO: this is broken... fix it

u0 = rand;
u1 = rand;
h = rand;
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

savg = @(lam) (1 - lam)*s0 + lam*s1;
stheta = @(lam) (1 - theta)*s + theta*savg(lam);

u = @(lam) (1 - lam)*u0 + lam*u1;
p = @(lam) p0 + lam*(p1 - p0);
q = @(lam) p(lam)'*p(lam);
l = @(lam) sqrt(q(lam));

F0 = @(lam) u(lam) + h*stheta(lam)*l(lam);
Lams = -0.5:0.01:1.5;
F0s = arrayfun(F0, Lams);

[F0opt, lamopt, lam1, lam2, degenerate] = F0_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1, false);

fprintf('F0 = %0.16g\n', F0opt);

figure;
hold on;
plot(Lams, F0s);
if ~degenerate
    plot(lam1, F0(lam1), '+');
    plot(lam2, F0(lam2), '+');
end
plot(lamopt, F0(lamopt), '*');
xlim([min(Lams) max(Lams)]);
