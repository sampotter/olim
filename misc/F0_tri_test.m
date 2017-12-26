clear;

u0 = rand;
u1 = rand;
h = rand;
s = rand;
s0 = rand;
s1 = rand;
theta = 0.5;

fprintf(['u0 = %g, u1 = %g, h = %g, s = %g, s0 = %g, s1 = %g, theta ' ...
         '= %g\n'], u0, u1, h, s, s0, s1, theta);

p0 = randn(3, 1);
p1 = randn(3, 1);

fprintf('p0 = (%g, %g, %g)\n', p0(1), p0(2), p0(3));
fprintf('p1 = (%g, %g, %g)\n', p1(1), p1(2), p1(3));

du = u1 - u0;
dp = p1 - p0;
savg = (s0 + s1)/2;
stheta = (1 - theta)*s + theta*savg;
alpha = -du/(stheta*h);

u = @(lam) (1 - lam)*u0 + lam*u1;
p = @(lam) p0 + lam*dp;
q = @(lam) p(lam)'*p(lam);
l = @(lam) sqrt(q(lam));

F0 = @(lam) u(lam) + stheta*h*l(lam);

Lams = linspace(-0.5, 1.5, 151);
F0s = zeros(size(Lams), 'like', Lams);
for i = 1:length(F0s)
    F0s(i) = F0(Lams(i));
end

tmp = alpha^2 - norm(dp, 2)^2;
a = norm(dp, 2)^2*tmp;
b = dp'*p0*tmp;
c = alpha^2*norm(p0, 2)^2 - (dp'*p0)^2;

degenerate = false;
if b*b < a*c
    degenerate = true;
    if F0(0) <= F0(1)
        arglam = 0;
    else
        arglam = 1;
    end
    fprintf('degenerate problem\n');
    fprintf('arglam = %g\n', arglam);
else
    lhs = -b/a;
    rhs = sqrt(b^2 - a*c)/a;

    lam1 = lhs - rhs;
    lam2 = lhs + rhs;

    check = @(lam) abs(alpha*l(lam) - dp'*p(lam));
    tol = 1e-10;
    if check(lam1) < check(lam2)
        arglam = lam1;
    else
        arglam = lam2;
    end
    if arglam < 0 || 1 < arglam
        if F0(0) <= F0(1)
            arglam = 0;
        else
            arglam = 1;
        end
    end
    fprintf('check1 = %g, check2 = %g\n', check(lam1), check(lam2));
    fprintf('lam1 = %g, lam2 = %g, arglam = %g\n', lam1, lam2, arglam);
end

figure;
hold on;
plot(Lams, F0s);
if ~degenerate
    plot(lam1, F0(lam1), '+');
    plot(lam2, F0(lam2), '+');
end
plot(arglam, F0(arglam), '*');
xlim([min(Lams) max(Lams)]);