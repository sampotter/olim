function newton222mp1(u0, u1, u2, s0, s1, s2, h)
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
        s0 = 1;
    end
    if nargin < 5
        s1 = 1;
    end
    if nargin < 6
        s2 = 1;
    end
    if nargin < 7
        h = 1;
    end

    du = [u1; u2] - u0;
    ds = [s1; s2] - s0;
    
    s = @(x, y) (1 - x - y)*s0 + x*s1 + y*s2;

    M = [2 1; 1 2];
    e = [-1; -1];
    f = 1;

    Q = @(x, y) [x y 1]*[M e; e' f]*[x; y; 1];
    l = @(x, y) sqrt((1 - x).^2 + (1 - y).^2 + (x + y).^2);
    
    g = @(x, y) du + h*(l(x, y)*ds + (s(x, y)/l(x, y))*(M*[x; y] + e));
    H = @(x, y) h*(2*(ds - (s(x, y)/(2*Q(x, y)))*(M*[x; y] + e))*(M*[x; ...
                        y] + e)' + s(x, y)*M)/l(x, y);

    X(1, :) = [1/3; 1/3];
    
    F_ = @(x, y) (1 - x - y)*u0 + x*u1 + y*u2 + h*s(x, y)*l(x, y);
    Fvals(1) = F_(X(1, 1), X(1, 2));

    niter = 200;
    step = 1;
    scale = 0.9;
    stepsizes(1) = 1;
    for k = (1:niter) + 1
        x = X(k - 1, 1);
        y = X(k - 1, 2);
        p = H(x, y)\g(x, y);

        % Backtracking line search
        tmp = [x; y] - step*p;
        while F_(tmp(1), tmp(2)) > Fvals(k - 1)
            step = scale*step;
            tmp = [x; y] - step*p;
        end
        stepsizes(k) = step;
        X(k, :) = tmp;

        Fvals(k) = F_(X(k, 1), X(k, 2));
        E(k - 1) = norm(p, 'inf')/norm(X(k, :), 'inf');
        if E(k - 1) < eps || abs(Fvals(k) - Fvals(k - 1)) < eps ...
                || any(X(k, :) < 0) || sum(X(k, :)) > 1
            break;
        end
    end
    if E(length(E)) == 0
        E(length(E)) = eps;
    end
    
    min_ = -1;
    max_ = 2;
    lin = linspace(min_, max_, 121);
    [x y] = meshgrid(lin, lin);
    F = (1 - x - y)*u0 + x*u1 + y*u2 + h*s(x, y).*l(x, y);

    figure;

    subplot(2, 2, 1);
    hold on;
    contour(x, y, F, 15);
    plot([0 1 0 0], [0 0 1 0], 'k');
    plot(X(:, 1), X(:, 2), '-*');
    xlim([min_, max_]);
    ylim([min_, max_]);
    scatter(X(size(X, 1), 1), X(size(X, 1), 2), 100);

    subplot(2, 2, 2);
    semilogy(1:length(E), E, '-x');
    xlim([0, length(E) + 1]);

    subplot(2, 2, 3);
    semilogy(1:length(Fvals), abs(Fvals - Fvals(length(Fvals))));
    xlim([1 length(Fvals)]);

    subplot(2, 2, 4);
    plot((1:length(stepsizes)) + 1, stepsizes);

end
