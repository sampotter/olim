function [U, lam] = tri_newton(u0, u1, s, s0, s1, h, method)
    if nargin < 1, u0 = 0.1; end
    if nargin < 2, u1 = 0; end
    if nargin < 3, s = 1; end
    if nargin < 4, s0 = 1.2; end
    if nargin < 5, s1 = 1.1; end
    if nargin < 6, h = 0.9; end
    if nargin < 7, method = '11'; end

    du = u1 - u0;
    sbar0 = (s + s0)/2;
    sbar1 = (s + s1)/2;
    dsbar = sbar1 - sbar0;

    if strcmp(method, '11')
        a = 2;
        b = -1;
        c = 1;
    elseif strcmp(method, '12')
        a = 1;
        b = 0;
        c = 1;
    elseif strcmp(method, '13')
        a = 2;
        b = 0;
        c = 1;
    elseif strcmp(method, '22')
        a = 2;
        b = -1;
        c = 2;
    elseif strcmp(method, '23')
        a = 1;
        b = 0;
        c = 2;
    end

    % function definitions
    u = @(x) (1 - x)*u0 + x*u1;
    sbar = @(x) (1 - x)*sbar0 + x*sbar1;
    q = @(x) a*x.^2 + 2*b*x + c;
    dq = @(x) 2*(a*x + b);
    l = @(x) sqrt(q(x));
    dl = @(x) dq(x)./(2*l(x));
    d2l = @(x) (a*c - b*b)./(q(x).*l(x));
    f = @(x) u(x) + h*sbar(x).*l(x);
    df = @(x) du + h*(dsbar*l(x) + sbar(x).*dl(x));
    d2f = @(x) h*(2*dsbar*dl(x) + sbar(x).*d2l(x));
    % p = @(x) -df(x)./d2f(x);
    p = @(x) -df(x)./max(0.1, d2f(x));

    x0 = 0.5;

    X(1) = x0;
    X(2) = X(1) + p(X(1));
    F(1) = f(X(1));
    F(2) = f(X(2));
    k = 2;
    while abs(F(k) - F(k - 1)) > eps && ...
            abs(p(X(k)))/abs(X(k)) > eps && ...
            0 <= X(k) && X(k) <= 1
        k = k + 1;
        X(k) = X(k - 1) + p(X(k - 1));
        F(k) = f(X(k));
    end
    if X(k) < 0 || 1 < X(k)
        lam = nan;
        U = inf;
    else
        lam = X(length(X));
        U = F(length(F));
        fprintf('df* = %0.16g\n', df(lam));
        fprintf('d2f* = %0.16g\n', d2f(lam));
    end
    
    K = length(F);
    assert(K == length(X));

    figure;
    subplot(2, 2, 1);
    xs = linspace(-0.5, 1.5, 201);
    plot(xs, f(xs)); hold on;
    plot(xs, df(xs));
    plot(xs, d2f(xs));
    plot([-0.5 1.5], [0 0], 'k--');
    legend('f', 'df', 'd2f');
    xlim([-0.5 1.5]);
    scatter(X, F, [], linspace(0, 1, K));
    subplot(2, 2, 2);
    plot(1:K, X, '-*');
    ylim([-0.5 1.5]);
    xlim([0 K+1])
    subplot(2, 2, 3);
    semilogy(2:K, abs(F(2:K) - F(1:K-1)), '-o');
    ylim([1e-18 1])
    xlim([1 K+1])
    subplot(2, 2, 4);
    semilogy(2:K, abs(X(2:K) - X(1:K-1)), '-o');
    ylim([1e-18 1])
    xlim([1 K+1])
    
    function alpha = get_alpha(k)
        alpha = 1;
        scale = 0.9;
        c1 = 1e-4;
        c2 = 1e-3;
        while f(X(k - 1) + alpha*p(X(k - 1))) > f(X(k - 1)) + ...
                c1*alpha*df(X(k - 1))*p(X(k - 1)) & ...
                abs(df(X(k - 1) + alpha*p(X(k - 1))*p(X(k - 1)))) > ...
                c2*abs(df(X(k - 1))*p(X(k - 1)))
            alpha = scale*alpha;
        end
    end
end
