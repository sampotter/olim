clear;

path(path, '../build/Release');

marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};

s1 = @(x, y, z) 1 - sin(sqrt(x.^2 + y.^2 + z.^2));
f1 = @(x, y, z) cos(sqrt(x.^2 + y.^2 + z.^2)) + sqrt(x.^2 + y.^2 + z.^2) - 1;

S = {s1};
F = {f1};

minMpower = 2;
maxMpower = 6;
Ms = (2.^(minMpower:maxMpower)) + 1;

n = 1;
s = S{n};
f = F{n};

methodnames = {'basic', 'olim6_rhr_arma', 'olim18_rhr_arma', 'olim26_rhr_arma'};

U = {};
E_inf = {};
E_2 = {};

relerr = @(method, u, U, p) norm(u(:) - U(:), p)/norm(u(:), p);

k = 1;
for M = Ms
    fprintf('M = %d\n', M);
    
    B = zeros(M, M, M, 'logical');
    B((M + 1)/2, (M + 1)/2, (M + 1)/2) = 1;
    h = 2/(M - 1);
    [X Y Z] = meshgrid(linspace(-1, 1, M), linspace(-1, 1, M), linspace(-1, 1, M));

    get_U = @(method) fmm(B, 'h', h, 'Speed', s, 'Method', method, ...
                          'x0', 1, 'y0', 1, 'z0', 1);

    u = f(X, Y, Z);
    u(isnan(u)) = 0;
    
    for k = 1:length(methodnames)
        method = methodnames{k};
        U{k} = get_U(method);
        E_inf{k} = relerr(method, u, U{k}, 'inf');
        E_2{k} = relerr(method, u, U{k}, 2);
    end
    
    k = k + 1;
end

figure;
set(gcf, 'Name', 'Relative Error', 'NumberTitle', 'off');

getplotsymb = @(index) strcat('-', marks{index});

subplot(1, 2, 1);
for k = 1:length(methodnames)
    method = methodnames{k}
    loglog(Ms, getfield(E_inf, method), getplotsymb(k));
    hold on;
end
title('l_\infty');
legend(methodnames);

subplot(1, 2, 2);
for k = 1:length(methodnames)
    method = methodnames{k}
    loglog(Ms, getfield(E_2, method), getplotsymb(k));
    hold on;
end
title('l_2');
legend(methodnames);