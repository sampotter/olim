clear;

path(path, '../build/Release');

s1 = @(x, y) 1 - sin(sqrt(x.^2 + y.^2));
f1 = @(x, y) cos(sqrt(x.^2 + y.^2)) + sqrt(x.^2 + y.^2) - 1;

s2 = @(x, y) abs(x + y);
f2 = @(x, y) ((x + y).^2)/(2*sqrt(2));

s3 = @(x, y) sqrt(x.^2 + y.^2);
f3 = @(x, y) (x.^2 + y.^2)/2;

% Breaks olim8_mp0 and mp1?
s4 = @(x, y) 2*sqrt((x.^2 - y.^2).^2.*(x.^4 + 14.*x.^2.*y.^2 + y.^4)./((x.^2 + y.^2).^3));
f4 = @(x, y) -3*x.^2 + y.^2 + 4*x.^4./(x.^2 + y.^2);

s5 = @(x, y) sqrt(x.^18 + y.^18);
f5 = @(x, y) (x.^10 + y.^10)/10;

s6 = @(x, y) 4*sqrt((x - y).^2.*(x + y).^2.*(x.^2 + y.^2));
f6 = @(x, y) (x.^2 - y.^2).^2;

s = s2;
f = f2;

k = 1;
Ms = 2*ceil(logspace(1, 3, 10)/2) + 1;
for M = Ms
    fprintf('M = %d\n', M);
    B = zeros(M, 'logical');
    B((M + 1)/2, (M + 1)/2) = 1;
    h = 2/(M - 1);
    [X Y] = meshgrid(linspace(-1, 1, M), linspace(-1, 1, M));
    u = f(X, Y);
    u(isnan(u)) = 0;
    U_basic = fmm(B, 'h', h, 'Speed', s, 'Method', 'basic', 'x0', 1, 'y0', 1);
    U_olim8_rhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);
    U_olim8_mp0 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp0', 'x0', 1, 'y0', 1);
    U_olim8_mp1 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp1', 'x0', 1, 'y0', 1);
    relerr = @(U) norm(u(:) - U(:), 'inf')/norm(u(:), 'inf');
    E_basic(k) = relerr(U_basic);
    E_olim8_rhr(k) = relerr(U_olim8_rhr);
    E_olim8_mp0(k) = relerr(U_olim8_mp0);
    E_olim8_mp1(k) = relerr(U_olim8_mp1);
    k = k + 1;
end

figure;
loglog(Ms, E_basic); hold on;
loglog(Ms, E_olim8_rhr); hold on;
loglog(Ms, E_olim8_mp0); hold on;
loglog(Ms, E_olim8_mp1); hold on;
xlim([min(Ms), max(Ms)]);
legend('basic', 'olim8\_rhr', 'olim8\_mp0', 'olim8\_mp1');
