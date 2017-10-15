function newton222rhr(u0, u1, u2, s, h);
    if nargin < 1
        u0 = 0;
    end
    if nargin < 2
        u1 = 0;
    end
    if nargin < 3
        u2 = 0;
    end
    if nargin < 4
        s = 1;
    end
    if nargin < 5
        h = 1;
    end
    
    du = [u1; u2] - u0;
    sh = s*h;

    l = @(x, y) sqrt((1 - x).^2 + (1 - y).^2 + (x + y).^2);
    FHessInv = @(x, y) (l(x, y)/sh)*[3 + x.*(3*x - 2), x.*(3*y - 1) - y ...
                        - 1; x.*(3*y - 1) - y - 1, 3 + y.*(3*y - 2)]/4;
    FGrad = @(x, y) du + (sh/l(x, y))*[2*x + y - 1; x + 2*y - 1];

    X(1, :) = [1/3; 1/3];

    for k = 2:50
        x = X(k - 1, 1);
        y = X(k - 1, 2);
        
        step = FHessInv(x, y)*FGrad(x, y);
        alpha = get_step_size([x; y], step, du, sh);
        X(k, :) = [x; y] + alpha*step;

        relerr = norm(X(k, :) - X(k - 1, :), 'inf')/norm(X(k, :));
        E(k - 1) = relerr;
        if relerr < 1e-15
            break;
        end
    end
    if E(length(E)) == 0
        E(length(E)) = eps; % for better plotting
    end

    lin = linspace(-2, 3, 251);
    [x_ y_] = meshgrid(lin, lin);
    F = (1 - x_ - y_)*u0 + x_*u1 + y_*u2 + s*h*l(x_, y_);

    figure;
    subplot(1, 2, 1);
    hold on;
    contour(x_, y_, F, 10);
    plot([0 1 0 0], [0 0 1 0], 'k');
    plot(X(:, 1), X(:, 2), '-*');
    scatter(X(size(X, 1), 1), X(size(X, 1), 2), 100);
    subplot(1, 2, 2);
    semilogy(1:length(E), E);
end

function e = get_e()
    e = [-2; -2];
end

function alpha = get_step_size(lambda, p, du, sh)
    e = get_e();
    delta = (2*du'*p/sh)^2;
    tmp1 = 2*Mdot(p, p) - delta;
    tmp2 = Mdot(lambda, p) + e'*p;
    a = Mdot(p, p)*tmp1/2
    b = tmp1*tmp2
    c = tmp2*tmp2 - delta*Q(lambda)
    disc = b*b - 4*a*c;
    assert(disc >= 0);
    lhs = -b/(2*a);
    rhs = sqrt(disc)/(2*a);
    if lhs + rhs <= 0
        alpha = lhs + rhs;
    else
        alpha = lhs - rhs;
    end
end

function value = Mdot(u, v)
    M = [4 2; 2 4];
    value = u'*M*v;
end

function value = Q(lambda)
    x = lambda(1);
    y = lambda(2);
    value = (1 - x).^2 + (1 - y).^2 + (x + y).^2;
end
