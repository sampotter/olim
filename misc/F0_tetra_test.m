clear;

F0_tetra_setup;

F0 = @(x) u(x) + sh*l(x);
dF0 = @(x) du + sh*dP'*p(x)/l(x);
d2F0 = @(x) sh*dP'*cprojp(x)*dP/l(x);

beta = @(x) ((pinv(dP)'*du)'*cprojp(x)*(pinv(dP)'*du))/(sh*sh);
stepsize = @(x) -1/sqrt(1 - beta(x));

out = sqp(F0, dF0, d2F0, -A, -b, ones(N - 1, 1)/N, eps, 100);

[X Y] = meshgrid(linspace(-0.2, 1.2, 141), linspace(-0.2, 1.2, 141));
Z = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z(i, j) = F0([X(i, j) Y(i, j)]);
    end
end

ndiv = 5;
W = plot_winding_number(dF0, ndiv);
fprintf('W = %g\n', W);

F0_plot;