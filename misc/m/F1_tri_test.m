clear;

u0 = 0.59457611725233628;
u1 = 0.56102682919361013;
du = u1 - u0;
s = 0.0025050133959455545;
s0 = 0.0025050133959455545;
s1 = 0.2382400185837108;
ds = s1 - s0;
h = 0.5;
theta = 0.5;
p0 = [1; 1; 0];
p1 = [0; 1; 1];
dp = p1 - p0;

% h = rand;
% u0 = rand;
% u1 = rand;
% s = rand;
% s0 = rand;
% s1 = rand;
% theta = rand;
% p0 = randn(3, 1);
% p1 = randn(3, 1);

fprintf(['u0 = %g, u1 = %g, h = %g, s = %g, s0 = %g, s1 = %g, theta ' ...
         '= %g\n'], u0, u1, h, s, s0, s1, theta);

u = @(lam) (1 - lam)*u0 + lam*u1;
s = @(lam) (1 - theta)*s + theta*((1 - lam)*s0 + lam*s1);
p = @(lam) p0 + lam*(p1 - p0);
l = @(lam) sqrt(p(lam)'*p(lam));
F1 = @(lam) u(lam) + h*s(lam)*l(lam);
dF1 = @(lam) du + h*(ds*l(lam)/2 + s(lam)*dp'*p(lam)/l(lam));

[lamopt, ~, deg] = hybrid(dF1, 0, 1, eps);
if deg
    if u0 + s(0)*h*l(0) < u1 + s(1)*h*l(1)
        lamopt = 0;
    else
        lamopt = 1;
    end
end
F1opt = F1(lamopt);

% [F1opt, lamopt, lam, F1iters] = F1_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1);

fprintf('F1 = %0.16g\nlamopt = %0.16g\n', F1opt, lamopt);

Lams = linspace(-0.5, 1.5, 151);
F1s = arrayfun(F1, Lams);

figure;
hold on;
plot(Lams, F1s);
% scatter(lam, F1iters, 'r');
plot(lamopt, F1(lamopt), '*');
xlim([min(Lams) max(Lams)]);
