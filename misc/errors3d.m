clear;

path(path, '../build/Release');

marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};

r = @(x, y, z) sqrt(x.^2 + y.^2 + z.^2);

s1 = true;
f1 = @(x, y, z) r(x, y, z);

s2 = @(x, y, z) 1 - sin(r(x, y, z));
f2 = @(x, y, z) cos(r(x, y, z)) + r(x, y, z) - 1;

S = {s1, s2};
F = {f1, f2};

n = 2;
s = S{n};
f = F{n};

ns = 2.^(2:6) + 1;
% ns = 5:2:31;

relerr = @(x, y, p) max(norm(x(:) - y(:), p)/norm(x(:), p), ...
                        norm(x(:) - y(:), p)/norm(y(:), p));

k = 1;
for n = ns
    fprintf('n = %d\n', n);

    % Compute grid
    B = zeros(n, n, n, 'logical');
    i0 = (n + 1)/2;
    B(i0, i0, i0) = 1;

    % Compute other parameters
    h = 2/(n - 1);
    
    % Compute ground truth solution
    L = linspace(-1, 1, n);
    [x y z] = meshgrid(L, L, L);
    u = f(x, y, z);
    
    % Compute solutions using different methods
    if islogical(s) && s
        Ubasic = fmm(B, 'h', h, 'Method', 'basic', 'x0', 1, 'y0', 1, 'z0', 1);
        U6mp0 = fmm(B, 'h', h, 'Method', 'olim6_mp0', 'x0', 1, 'y0', 1, 'z0', 1);
        U6rhr = fmm(B, 'h', h, 'Method', 'olim6_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
        U18mp0 = fmm(B, 'h', h, 'Method', 'olim18_mp0', 'x0', 1, 'y0', 1, 'z0', 1);
        U18rhr = fmm(B, 'h', h, 'Method', 'olim18_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
        U26 = fmm(B, 'h', h, 'Method', 'olim26_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
    else
        Ubasic = fmm(B, 'h', h, 'Speed', s, 'Method', 'basic', 'x0', 1, 'y0', 1, 'z0', 1);
        U6mp0 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim6_mp0', 'x0', 1, 'y0', 1, 'z0', 1);
        U6rhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim6_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
        U18mp0 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim18_mp0', 'x0', 1, 'y0', 1, 'z0', 1);
        U18rhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim18_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
        U26 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim26_rhr', 'x0', 1, 'y0', 1, 'z0', 1);
    end

    % Compute relative errors
    Ebasic(k) = relerr(Ubasic, u, 'inf');
    E6mp0(k) = relerr(U6mp0, u, 'inf');
    E6rhr(k) = relerr(U6rhr, u, 'inf');
    E18mp0(k) = relerr(U18mp0, u, 'inf');
    E18rhr(k) = relerr(U18rhr, u, 'inf');
    E26(k) = relerr(U26, u, 'inf');
    
    k = k + 1;
end

getplotsymb = @(index) strcat('-', marks{index});

% Plot errors
figure;
loglog(ns, Ebasic, getplotsymb(1));
hold on;
loglog(ns, E6mp0, getplotsymb(2));
loglog(ns, E6rhr, getplotsymb(3));
loglog(ns, E18mp0, getplotsymb(4));
loglog(ns, E18rhr, getplotsymb(5));
loglog(ns, E26, getplotsymb(6));
ylabel('||u - U||_\inf/||u||_\inf');
xlabel('n (s.t. # nodes = n^3)');
legend('Basic', 'OLIM6 (mp0)', 'OLIM6 (rhr)', 'OLIM18 (mp0)', ...
       'OLIM18 (rhr)', 'OLIM26 (rhr)');
