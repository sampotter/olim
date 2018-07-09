clear;

%%% TODO: this is broken... fix it

% u0 = rand;
% u1 = rand;
% du = u1 - u0;
% h = 1;
% s = rand;
% s0 = rand;
% s1 = rand;
% theta = rand;
% p0 = rand(3, 1);
% p1 = rand(3, 1);
% dp = p1 - p0;
% rfac = 0.1;
% pfac = (rfac/h)*randn(3, 1);
% sfac = rand;

u0 = 0.058596292634851235;
u1 = 0.058596292634851235;
du = u1 - u0;
h = 0.0625;
s = 0.9117266963977384;
s0 = 0.93754068215761976;
s1 = 0.93754068215761976;
p0 = [1; 0];
p1 = [0; -1];
dp = p1 - p0;
pfac = [1; -1];
sfac = 1;

% fprintf(['u0 = %g, u1 = %g, h = %g, s = %g, s0 = %g, s1 = %g, theta ' ...
%          '= %g\n'], u0, u1, h, s, s0, s1, theta);
% fprintf('p0 = (%g, %g, %g)\n', p0(1), p0(2), p0(3));
% fprintf('p1 = (%g, %g, %g)\n', p1(1), p1(2), p1(3));

u = @(lam) (1 - lam)*u0 + lam*u1;
p = @(lam) p0 + lam*(p1 - p0);
q = @(lam) p(lam)'*p(lam);
l = @(lam) sqrt(q(lam));
n = @(lam) p(lam)/l(lam);

lfac = @(lam) norm(p(lam) - pfac);
nfac = @(lam) (p(lam) - pfac)/lfac(lam);
T = @(lam) sfac*h*lfac(lam);

tau0 = u0 - T(0);
tau1 = u1 - T(1);
dtau = tau1 - tau0;

tau = @(lam) tau0 + dtau*lam;

F0 = @(lam) u(lam) + s*h*l(lam);
dF0 = @(lam) du + s*h*dot(dp, n(lam));

F0f = @(lam) tau(lam) + T(lam) + s*h*l(lam);
dF0f = @(lam) dtau + sfac*h*dot(dp, nfac(lam)) + s*h*dot(dp, n(lam));

% [F0opt, lamopt, lam1, lam2, degenerate] = F0_tri_newton(u0, u1,
% h, s, s0, s1, theta, p0, p1, false);

[lamopt, ~, deg] = hybrid(dF0, 0, 1, eps);
if deg
    if F0(0) < F0(1)
        lamopt = 0;
        F0opt = F0(0);
    else
        lamopt = 1;
        F0opt = F0(1);
    end
end
F0opt = F0(lamopt);
fprintf('lamopt = %0.16g\n', lamopt);
fprintf('F0 = %0.16g\n', F0opt);

[lamoptf, ~, deg] = hybrid(dF0f, 0, 1, eps);
if deg
    if F0f(0) < F0f(1)
        lamoptf = 0;
        F0optf = F0f(0);
    else
        lamoptf = 1;
        F0optf = F0f(1);
    end
end
F0optf = F0f(lamoptf);
fprintf('lamoptf = %0.16g\n', lamoptf);
fprintf('F0f = %0.16g\n', F0optf);

figure;

subplot(121);
hold on;
Lams = -0.5:0.01:1.5;
plot(Lams, arrayfun(F0, Lams), '-r');
plot(Lams, arrayfun(F0f, Lams), '-b');
plot(lamopt, F0(lamopt), '*r');
plot(lamoptf, F0f(lamoptf), '*b');
xlim([min(Lams) max(Lams)]);

subplot(122);
hold on;
Lams = -0.5:0.01:1.5;
plot(Lams, arrayfun(@(x) F0(x) - F0f(x), Lams));
xlim([min(Lams) max(Lams)]);