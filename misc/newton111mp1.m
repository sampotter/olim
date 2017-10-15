function newton111mp1(u0, u1, u2, s0, s1, s2, h)
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
    ds = [s1; s2] - s0;
    
    s = @(x, y) (1 - x - y)*s0 + x*s1 + y*s2;

    M = [2 1; 1 2];
    e = [-1; -1];
    f = 1;

    l = @(x, y) sqrt((1 - x - y).^2 + x.^2 + y.^2);
    Q = @(x, y) [x y 1]*[M e; e' f]*[x; y; 1];
    
    g = @(x, y) du + h*(l(x, y)*ds + (s(x, y)/l(x, y))*(M*[x; y] + e));
    H = @(x, y) h*(2*(ds - (s(x, y)/(2*Q(x, y)))*(M*[x; y] + e))*(M*[x; ...
                        y] + e)' + s(x, y)*M)/l(x, y);

    X(1, :) = [1/3; 1/3];
    
    niter = 50;
    for k = (1:niter) + 1
        x = X(k - 1, 1);
        y = X(k - 1, 2);
        p = H(x, y)\g(x, y);
        X(k, :) = [x; y] - p;
        E(k - 1) = norm(p, 'inf')/norm(X(k, :), 'inf');
        if E(k - 1) < 1e-15
            break;
        end
    end
    if E(length(E)) == 0
        E(length(E)) = eps;
    end
    
    lin = linspace(-0.1, 1.1, 121);
    [x y] = meshgrid(lin, lin);
    F = (1 - x - y)*u0 + x*u1 + y*u2 + h*s(x, y).*l(x, y);

    figure;
    subplot(1, 2, 1);
    hold on;
    contour(x, y, F, 15);
    plot([0 1 0 0], [0 0 1 0], 'k');
    plot(X(:, 1), X(:, 2), '-*');
    scatter(X(size(X, 1), 1), X(size(X, 1), 2), 100);
    subplot(1, 2, 2);
    semilogy(1:length(E), E, '-x');
    xlim([0, length(E) + 1]);
end
