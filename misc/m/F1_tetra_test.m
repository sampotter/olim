clear;

F0_tetra_setup;

stheta = @(x) (1 - theta)*s + theta*((1 - sum(x))*s0 + x(1)*s1 + x(2)*s2);

F1 = @(x) u(x) + h*stheta(x)*l(x);
dF1 = @(x) du + h*(theta*q(x)*ds + stheta(x)*dP'*p(x))/l(x);
d2F1 = @(x) h*(theta*(dP'*p(x)*ds' + ds*p(x)'*dP) + stheta(x)*dP'*cprojp(x)*dP)/l(x);

out = sqp(F1, dF1, d2F1, -A, -b, ones(N - 1, 1)/N, eps, 100);
fprintf('U = %g\n', out.fs(out.iters));
lamopt = out.xs(:, out.iters + 1);

[X Y] = meshgrid(linspace(-0.2, 1.2, 141), linspace(-0.2, 1.2, 141));
Z = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z(i, j) = F1([X(i, j) Y(i, j)]);
    end
end

[~, lam01, ~, ~] = F1_tri_newton(u0, u1, h, s, s0, s1, theta, p0, p1);
[~, lam02, ~, ~] = F1_tri_newton(u0, u2, h, s, s0, s2, theta, p0, p2);
[~, lam12, ~, ~] = F1_tri_newton(u1, u2, h, s, s1, s2, theta, p1, p2);

lam01 = [lam01; 0];
lam02 = [0; lam02];
lam12 = [1 - lam12; lam12];

mu = @(x, I) (A(I, :)*inv(d2F1(x))*A(I, :)')\(A(I, :)*inv(d2F1(x))*dF1(x));

I01 = find(abs(A*lam01 - b) < eps);
I02 = find(abs(A*lam02 - b) < eps);
I12 = find(abs(A*lam12 - b) < eps);

mu01 = mu(lam01, I01);
mu02 = mu(lam02, I02);
mu12 = mu(lam12, I12);

for i = 1:length(mu01)
    fprintf('mu01(%d) = %g\n', i, mu01(i));
end
for i = 1:length(mu02)
    fprintf('mu02(%d) = %g\n', i, mu02(i));
end
for i = 1:length(mu12)
    fprintf('mu12(%d) = %g\n', i, mu12(i));
end

% ndiv = 5;
% W = plot_winding_number(dF1, ndiv);
% fprintf('W = %g\n', W);

plot_tri_optima = true;
using_Z0_and_Z1 = false;

F0_plot;
