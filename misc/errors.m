clear;

path(path, '../build/Release');

s1 = @(x, y) 1 - sin(sqrt(x.^2 + y.^2));
f1 = @(x, y) cos(sqrt(x.^2 + y.^2)) + sqrt(x.^2 + y.^2) - 1;

s2 = @(x, y) abs(x + y);
f2 = @(x, y) ((x + y).^2)/(2*sqrt(2));

s3 = @(x, y) sqrt(x.^2 + y.^2);
f3 = @(x, y) (x.^2 + y.^2)/2;

% Breaks olim8_mp0 and mp1? (because of singularity at origin)
s4 = @(x, y) 2*sqrt((x.^2 - y.^2).^2.*(x.^4 + 14.*x.^2.*y.^2 + y.^4)./((x.^2 + y.^2).^3));
f4 = @(x, y) -3*x.^2 + y.^2 + 4*x.^4./(x.^2 + y.^2);

s5 = @(x, y) sqrt(x.^18 + y.^18);
f5 = @(x, y) (x.^10 + y.^10)/10;

s6 = @(x, y) 4*sqrt((x - y).^2.*(x + y).^2.*(x.^2 + y.^2));
f6 = @(x, y) (x.^2 - y.^2).^2;

s7 = @(x, y) masha_s(x, y);
f7 = @(x, y) masha_f(x, y);

s8 = @(x, y) (1 + sqrt(x.^2 + y.^2)).^-2;
f8 = @(x, y) 1 - 1./(1 + sqrt(x.^2 + y.^2));

a = 10;
s9 = @(x, y) 2.*a.*sqrt(exp(-2.*a.*(x.^2 + y.^2)).*(x.^2 + y.^2));
f9 = @(x, y) 1 - 1./exp(a.*(x.^2 + y.^2));

S = {s1 s2 s3 s4 s5 s6 s7 s8 s9};
F = {f1 f2 f3 f4 f5 f6 f7 f8 f9};

n = 1;
s = S{n};
f = F{n};

k = 1;
% Ms = 2*ceil(logspace(1, 3, 11)/2) + 1;
Ms = (2.^(3:7)) + 1;
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
    % u = sqrt(X.*X + Y.*Y);
    % U_basic = fmm(B, 'h', h, 'Method', 'basic', 'x0', 1, 'y0', 1);
    % U_olim8_rhr = fmm(B, 'h', h, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);
    % U_olim8_mp0 = fmm(B, 'h', h, 'Method', 'olim8_mp0', 'x0', 1, 'y0', 1);
    % U_olim8_mp1 = fmm(B, 'h', h, 'Method', 'olim8_mp1', 'x0', 1, 'y0', 1);
    relerr = @(U, p) norm(u(:) - U(:), p)/norm(u(:), p);
    E_basic_inf(k) = relerr(U_basic, 'inf');
    E_olim8_rhr_inf(k) = relerr(U_olim8_rhr, 'inf');
    E_olim8_mp0_inf(k) = relerr(U_olim8_mp0, 'inf');
    E_olim8_mp1_inf(k) = relerr(U_olim8_mp1, 'inf');
    E_basic_2(k) = relerr(U_basic, 2);
    E_olim8_rhr_2(k) = relerr(U_olim8_rhr, 2);
    E_olim8_mp0_2(k) = relerr(U_olim8_mp0, 2);
    E_olim8_mp1_2(k) = relerr(U_olim8_mp1, 2);
    k = k + 1;
end

figure;
set(gcf, 'Name', 'Relative Error', 'NumberTitle', 'off');

subplot(1, 2, 1);
loglog(Ms, E_basic_inf); hold on;
loglog(Ms, E_olim8_rhr_inf); hold on;
loglog(Ms, E_olim8_mp0_inf); hold on;
loglog(Ms, E_olim8_mp1_inf); hold on;
title('l_\infty');
ylabel('||u - U||_\infty/||u||_\infty');
xlabel('n');
xlim([min(Ms), max(Ms)]);
legend('basic', 'olim8\_rhr', 'olim8\_mp0', 'olim8\_mp1');
subplot(1, 2, 2);
loglog(Ms, E_basic_2); hold on;
loglog(Ms, E_olim8_rhr_2); hold on;
loglog(Ms, E_olim8_mp0_2); hold on;
loglog(Ms, E_olim8_mp1_2); hold on;
title('l_2');
ylabel('||u - U||_2/||u||_2');
xlabel('n');
xlim([min(Ms), max(Ms)]);
legend('basic', 'olim8\_rhr', 'olim8\_mp0', 'olim8\_mp1');

figure;
set(gcf, 'Name', 'Analytic Solution', 'NumberTitle', 'off');

imagesc(u);
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

figure;
set(gcf, 'Name', 'Pointwise Error', 'NumberTitle', 'off');

subplot(2, 2, 1); 
imagesc(U_basic - u);
title('basic'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(2, 2, 2); 
imagesc(U_olim8_rhr - u);
title('olim8\_rhr'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(2, 2, 3); 
imagesc(U_olim8_mp0 - u);
title('olim8\_mp0'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(2, 2, 4); 
imagesc(U_olim8_mp1 - u);
title('olim8\_mp1'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;
