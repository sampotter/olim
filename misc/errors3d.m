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

k = 1;
for M = Ms
    fprintf('M = %d\n', M);
    
    B = zeros(M, M, M, 'logical');
    B((M + 1)/2, (M + 1)/2, (M + 1)/2) = 1;
    h = 2/(M - 1);
    [X Y Z] = meshgrid(linspace(-1, 1, M), linspace(-1, 1, M), linspace(-1, 1, M));

    u = f(X, Y, Z);
    u(isnan(u)) = 0;
    
    U_olim6_rhr_arma = fmm(B, 'h', h, 'Speed', s, 'Method', ...
                           'olim6_rhr_arma', 'x0', 1, 'y0', 1, 'z0', 1);
    
    relerr = @(U, p) norm(u(:) - U(:), p)/norm(u(:), p);

    E_olim6_rhr_arma_inf(k) = relerr(U_olim6_rhr_arma, 'inf');
    
    E_olim6_rhr_arma_2(k) = relerr(U_olim6_rhr_arma, 2);
    
    k = k + 1;
end

figure;
set(gcf, 'Name', 'Relative Error', 'NumberTitle', 'off');

getplotsymb = @(index) strcat('-', marks{index});

subplot(1, 2, 1);
loglog(Ms, E_olim6_rhr_arma_inf, getplotsymb(1)); hold on;

subplot(1, 2, 2);
loglog(Ms, E_olim6_rhr_arma_2, getplotsymb(1)); hold on;