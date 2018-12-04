clear;

F0_tetra_setup;

F0 = @(x) u(x) + sh*l(x);
dF0 = @(x) du + sh*dP'*p(x)/l(x);
d2F0 = @(x) sh*dP'*cprojp(x)*dP/l(x);

F0fac = @(x) tau(x) + T(x) + sh*l(x);
dF0fac = @(x) dtau + sfac*h*dP'*nfac(x) + sh*dP'*n(x);
d2F0fac = @(x) d2F0(x) + ...
          (lfac(x) > 1e1*eps)*(sfac*h*dP'*(I - nfac(x)*nfac(x)')*dP/(max(eps, lfac(x))));

sh = @(x) h*((1 - theta)*s + theta*((1 - sum(x))*s0 + x(1)*s1 + x(2)*s2));
F1 = @(x) u(x) + sh(x)*l(x);
dF1 = @(x) du + (h*theta*q(x)*ds + sh(x)*dP'*p(x))/l(x);
d2F1 = @(x) (h*theta*(dP'*p(x)*ds' + ds*p(x)'*dP) + sh(x)*dP'*cprojp(x)*dP)/l(x);

beta = @(x) ((pinv(dP)'*du)'*cprojp(x)*(pinv(dP)'*du))/(sh*sh);
stepsize = @(x) -1/sqrt(1 - beta(x));

lam0 = ones(N - 1, 1)/N;

out0 = sqp(F0, dF0, d2F0, -A, -b, lam0, eps, 100);
out1 = sqp(F1, dF1, d2F1, -A, -b, lam0, eps, 100);

lam0gt = out0.xs(:, out0.iters + 1);
lam1gt = out1.xs(:, out1.iters + 1);

% fprintf('F1(lam0gt) = %0.16g\n', F1(lam0gt));

out0fac = sqp(F0fac, dF0fac, d2F0fac, -A, -b, lam0, eps, 100);

lam0facgt = out0fac.xs(:, out0fac.iters + 1);

% Try our unconstrained approach here...

lam0_in_bounds = abs(max(A*lam0gt - b)) > 10*eps;
lam1_in_bounds = abs(max(A*lam1gt - b)) > 10*eps;

% if lam0_in_bounds && lam1_in_bounds
%     fprintf('Trying unconstrained approach...\n')

%     [Q R] = qr(dP, 0);
%     sh_ = sh([1; 1]/3);
%     alpha = sqrt((p0'*(eye(3) - Q*Q')*p0)/(1 - norm(R'\(du/sh_))^2));
%     lam0 = -inv(R)*(alpha*inv(R')*(du/sh_) + Q'*p0);

%     k = 0;
%     lam(:, 1) = lam0;
%     g = @(x) -d2F1(x)\dF1(x);
%     while norm(g(lam(:, k + 1)), 'inf') > eps
%         k = k + 1;
%         lam(:, k + 1) = lam(:, k) + g(lam(:, k));
%     end
%     k = k + 1;
% end

% ---

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

Z0fac = zeros(size(X), 'like', X);
for i = 1:size(X)
    for j = 1:size(Y)
        Z0fac(i, j) = F0fac([X(i, j) Y(i, j)]);
    end
end

% ndiv = 5;
% W = plot_winding_number(dF0, ndiv);
% fprintf('W = %g\n', W);

using_Z0_and_Z1 = true;

F0_plot;

% Taking a look at skipping updates using KKT theory:

mu0 = @(x, I) (A(I, :)*inv(d2F0(x))*A(I, :)')\(A(I, :)*inv(d2F0(x))*dF0(x));
mu1 = @(x, I) (A(I, :)*inv(d2F1(x))*A(I, :)')\(A(I, :)*inv(d2F1(x))*dF1(x));

mu0p0 = mu0([0; 0], [1; 2]);
mu0p1 = mu0([1; 0], [2; 3]);
mu0p2 = mu0([0; 1], [1; 3]);

mu1p0 = mu1([0; 0], [1; 2]);
mu1p1 = mu1([1; 0], [2; 3]);
mu1p2 = mu1([0; 1], [1; 3]);

fprintf('mu0p0 = [%g; %g]\n', mu0p0);
fprintf('mu0p1 = [%g; %g]\n', mu0p1);
fprintf('mu0p2 = [%g; %g]\n', mu0p2);

fprintf('mu1p0 = [%g; %g]\n', mu1p0);
fprintf('mu1p1 = [%g; %g]\n', mu1p1);
fprintf('mu1p2 = [%g; %g]\n', mu1p2);

% Check solutions with predicted skips

F0_interior_soln = all(abs(A*lam0gt - b) > 10*eps);
F1_interior_soln = all(abs(A*lam1gt - b) > 10*eps);

skip_F0_tetra = any(mu0p0 < 0) || any(mu0p1 < 0) || any(mu0p2 < 0);
skip_F1_tetra = any(mu1p0 < 0) || any(mu1p1 < 0) || any(mu1p2 < 0);

fprintf('F0 skip = %d... OK? %d\n', skip_F0_tetra, ~F0_interior_soln);
fprintf('F1 skip = %d... OK? %d\n', skip_F1_tetra, ~F1_interior_soln);

