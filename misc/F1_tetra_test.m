clear;

F0_tetra_setup;

stheta = @(x) (1 - theta)*s + theta*((1 - sum(x))*s0 + x(1)*s1 + x(2)*s2);

F1 = @(x) u(x) + h*stheta(x)*l(x);
dF1 = @(x) du + h*(theta*q(x)*ds + stheta(x)*dP'*p(x))/l(x);
d2F1 = @(x) h*(theta*(dP'*p(x)*ds' + ds*p(x)'*dP) + stheta(x)*dP'*cprojp(x)*dP)/l(x);

out = sqp(F1, dF1, d2F1, -A, -b, ones(N - 1, 1)/N, eps, 100);

[X Y] = meshgrid(linspace(-0.2, 1.2, 141), linspace(-0.2, 1.2, 141));
Z = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z(i, j) = F1([X(i, j) Y(i, j)]);
    end
end

ndiv = 5;
W = plot_winding_number(dF1, ndiv);
fprintf('W = %g\n', W);

F0_plot;