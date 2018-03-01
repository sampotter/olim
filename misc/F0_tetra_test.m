clear;

F0_tetra_setup;

F0 = @(x) u(x) + sh*l(x);
dF0 = @(x) du + sh*dP'*p(x)/l(x);
d2F0 = @(x) sh*dP'*cprojp(x)*dP/l(x);

sh = @(x) h*((1 - theta)*s + theta*((1 - sum(x))*s0 + x(1)*s1 + x(2)*s2));
F1 = @(x) u(x) + sh(x)*l(x);
dF1 = @(x) du + (h*theta*q(x)*ds + sh(x)*dP'*p(x))/l(x);
d2F1 = @(x) (h*theta*(dP'*p(x)*ds' + ds*p(x)'*dP) + sh(x)*dP'*cprojp(x)*dP)/l(x);

beta = @(x) ((pinv(dP)'*du)'*cprojp(x)*(pinv(dP)'*du))/(sh*sh);
stepsize = @(x) -1/sqrt(1 - beta(x));

out0 = sqp(F0, dF0, d2F0, -A, -b, ones(N - 1, 1)/N, eps, 100);
out1 = sqp(F1, dF1, d2F1, -A, -b, ones(N - 1, 1)/N, eps, 100);

[X Y] = meshgrid(linspace(-0.2, 1.2, 141), linspace(-0.2, 1.2, 141));

Z0 = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z0(i, j) = F0([X(i, j) Y(i, j)]);
    end
end

Z1 = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z1(i, j) = F1([X(i, j) Y(i, j)]);
    end
end

% ndiv = 5;
% W = plot_winding_number(dF0, ndiv);
% fprintf('W = %g\n', W);

using_Z0_and_Z1 = true;

F0_plot;