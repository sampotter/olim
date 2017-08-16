clear;

path(path, '../build/Release');

marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};

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

s10 = @(x, y) star_s(x, y);

S = {s1 s2 s3 s4 s5 s6 s7 s8 s9, s10};
F = {f1 f2 f3 f4 f5 f6 f7 f8 f9, false};

n = 1;

maxMpower = 10;
Ms = (2.^(3:maxMpower)) + 1;

if n == 10
    fprintf('Loading speed function ground truth values...\n');
    if exist('U_s10_gt.mat', 'file') ~= 2
        M = pow2(nextpow2(maxMpower)) + 1;
        fprintf('Creating with M = %d...\n', M);
        make_s10_gt(M);
    end
    s10_gt = load('U_s10_gt.mat', 'M');
    M = s10_gt.M;
    if M < pow2(maxMpower) || round(log2(M - 1)) - log2(M - 1)
        make_s10_gt(pow2(nextpow2(maxMpower)) + 1);
    end
    s10_gt = load('U_s10_gt.mat');
end

s = S{n};
f = F{n};

k = 1;

for M = Ms
    fprintf('M = %d\n', M);

    B = zeros(M, 'logical');
    B((M + 1)/2, (M + 1)/2) = 1;
    h = 2/(M - 1);
    [X Y] = meshgrid(linspace(-1, 1, M), linspace(-1, 1, M));

    if n == 10
        step = round((s10_gt.M - 1)/(M - 1));
        I = 1:step:s10_gt.M;
        u = s10_gt.U(I, I);
    else
        u = f(X, Y);
    end
    u(isnan(u)) = 0;

    U_olim4_rhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim4_rhr', 'x0', 1, 'y0', 1);
    U_olim4_mp0 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim4_mp0', 'x0', 1, 'y0', 1);
    U_olim8_rhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);
    U_olim8_mp0 = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp0', 'x0', 1, 'y0', 1);
    U_olim8_mp1_bsearch = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp1_bsearch', ...
                      'x0', 1, 'y0', 1);
    U_olim8_mp1_gsl = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp1_gsl', ...
                      'x0', 1, 'y0', 1);

    % U_olim8_mp1_bsearch = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_mp1_bsearch', 'x0', 1, 'y0', 1);
    % u = sqrt(X.*X + Y.*Y);
    % U_olim8_rhr = fmm(B, 'h', h, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);
    % U_olim8_mp0 = fmm(B, 'h', h, 'Method', 'olim8_mp0', 'x0', 1, 'y0', 1);
    % U_olim8_mp1_bsearch = fmm(B, 'h', h, 'Method', 'olim8_mp1_bsearch', 'x0', 1,
    % 'y0', 1);

    relerr = @(U, p) norm(u(:) - U(:), p)/norm(u(:), p);

    E_olim4_rhr_inf(k) = relerr(U_olim4_rhr, 'inf');
    E_olim4_mp0_inf(k) = relerr(U_olim4_mp0, 'inf');
    E_olim8_rhr_inf(k) = relerr(U_olim8_rhr, 'inf');
    E_olim8_mp0_inf(k) = relerr(U_olim8_mp0, 'inf');
    E_olim8_mp1_bsearch_inf(k) = relerr(U_olim8_mp1_bsearch, 'inf');
    E_olim8_mp1_gsl_inf(k) = relerr(U_olim8_mp1_gsl, 'inf');

    E_olim4_rhr_2(k) = relerr(U_olim4_rhr, 2);
    E_olim4_mp0_2(k) = relerr(U_olim4_mp0, 2);
    E_olim8_rhr_2(k) = relerr(U_olim8_rhr, 2);
    E_olim8_mp0_2(k) = relerr(U_olim8_mp0, 2);
    E_olim8_mp1_bsearch_2(k) = relerr(U_olim8_mp1_bsearch, 2);
    E_olim8_mp1_gsl_2(k) = relerr(U_olim8_mp1_gsl, 2);

    k = k + 1;
end

figure;
set(gcf, 'Name', 'Relative Error', 'NumberTitle', 'off');

getplotsymb = @(index) strcat('-', marks{index});

subplot(1, 2, 1);
loglog(Ms, E_olim4_rhr_inf, getplotsymb(1)); hold on;
loglog(Ms, E_olim4_mp0_inf, getplotsymb(2)); hold on;
loglog(Ms, E_olim8_rhr_inf, getplotsymb(3)); hold on;
loglog(Ms, E_olim8_mp0_inf, getplotsymb(4)); hold on;
loglog(Ms, E_olim8_mp1_bsearch_inf, getplotsymb(5)); hold on;
loglog(Ms, E_olim8_mp1_gsl_inf, getplotsymb(6)); hold on;
title('l_\infty');
ylabel('||u - U||_\infty/||u||_\infty');
xlabel('n');
xlim([min(Ms), max(Ms)]);
legend('olim4\_rhr', ...
       'olim4\_mp0', ...
       'olim8\_rhr', ...
       'olim8\_mp0', ...
       'olim8\_mp1\_bsearch', ...
       'olim8\_mp1\_gsl');

subplot(1, 2, 2);
loglog(Ms, E_olim4_rhr_2, getplotsymb(1)); hold on;
loglog(Ms, E_olim4_mp0_2, getplotsymb(2)); hold on;
loglog(Ms, E_olim8_rhr_2, getplotsymb(3)); hold on;
loglog(Ms, E_olim8_mp0_2, getplotsymb(4)); hold on;
loglog(Ms, E_olim8_mp1_bsearch_2, getplotsymb(5)); hold on;
loglog(Ms, E_olim8_mp1_gsl_2, getplotsymb(6)); hold on;
title('l_2');
ylabel('||u - U||_2/||u||_2');
xlabel('n');
xlim([min(Ms), max(Ms)]);
legend('olim4\_rhr', ...
       'olim4\_mp0', ...
       'olim8\_rhr', ...
       'olim8\_mp0', ...
       'olim8\_mp1\_bsearch', ...
       'olim8\_mp1\_gsl');

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

subplot(3, 2, 1); 
imagesc(U_olim4_rhr - u);
title('olim4\_rhr'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(3, 2, 2); 
imagesc(U_olim4_mp0 - u);
title('olim4\_mp0'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(3, 2, 3); 
imagesc(U_olim8_rhr - u);
title('olim8\_rhr'); 
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(3, 2, 4); 
imagesc(U_olim8_mp0 - u);
title('olim8\_mp0');
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(3, 2, 5);
imagesc(U_olim8_mp1_bsearch - u);
title('olim8\_mp1\_bsearch');
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;

subplot(3, 2, 6);
imagesc(U_olim8_mp1_gsl - u);
title('olim8\_mp1\_gsl');
set(gca, 'XTick', [1 M/2 M]);
set(gca, 'XTickLabels', [-1 0 1]);
set(gca, 'YTick', [1 M/2 M]);
set(gca, 'YTickLabels', [-1 0 1]);
colorbar;
