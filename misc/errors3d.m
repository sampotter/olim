clear;

path(path, '../build/Release');

marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};

ns = 2.^(2:9) + 1;

relerr = @(x, y, p) max(norm(x(:) - y(:), p)/norm(x(:), p), ...
                        norm(x(:) - y(:), p)/norm(y(:), p));

k = 1;
for n = ns
    fprintf('n = %d\n', n);

    % Compute grid
    B = zeros(n, n, n, 'logical');
    B((n + 1)/2, (n + 1)/2, (n + 1)/2) = 1;

    % Compute other parameters
    h = 2/(n - 1);
    
    % Compute ground truth solution
    L = linspace(-1, 1, n);
    [x y z] = meshgrid(L, L, L);
    u = sqrt(x.^2 + y.^2 + z.^2);
    
    % Compute solutions using different methods
    Ubasic = fmm(B, 'h', h, 'Method', 'basic', 'x0', 1, 'y0', 1, 'z0', 1);
    U6 = fmm(B, 'h', h, 'Method', 'olim6_rhr_arma', 'x0', 1, 'y0', 1, 'z0', 1);

    Ebasic(k) = relerr(Ubasic, u, 'inf');
    E6(k) = relerr(U6, u, 'inf');
    
    k = k + 1;
end

figure;
loglog(ns, Ebasic, '*');
hold on;
loglog(ns, E6, '+');
legend('Basic', 'OLIM6 (right-hand rule)');