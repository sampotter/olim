clear;

path(path, '../build/Release');

B = zeros(21, 'logical');
B(11, 11) = 1;

h = 0.1;
s = @(x, y) masha_s(x, y);

Ubasic = fmm(B, 'h', h, 'Speed', s, 'Method', 'basic', 'x0', 1, 'y0', 1);

Ubasic = importdata('test2.txt');

ubasic = importdata('gfun_basic.txt');
Ibasic = ubasic < 1e5; % mask region with valid solution
ubasic = circshift(ubasic, -1, 2); % align these two so that
Ibasic = circshift(Ibasic, -1, 2); % the minimum is at the center NB:
                                   % the minimum in the file is != 0,
                                   % but more like 1e-7... okay for
                                   % testing, though


rebasic = NaN(21);
rebasic(Ibasic) = (Ubasic(Ibasic) - ubasic(Ibasic))./abs(ubasic(Ibasic));
rebasic(11, 11) = 0;

Urhr = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);

urhr = importdata('gfun_8pt_rhr.txt');
Irhr = urhr < 1e5; % mask region with valid solution
urhr = circshift(urhr, -1, 2); % align these two so that
Irhr = circshift(Irhr, -1, 2); % the minimum is at the center NB:
                                   % the minimum in the file is != 0,
                                   % but more like 1e-7... okay for
                                   % testing, though


rerhr = NaN(21);
rerhr(Irhr) = (Urhr(Irhr) - urhr(Irhr))./abs(urhr(Irhr));
rerhr(11, 11) = 0;

figure;
subplot(1, 2, 1);
imagesc(rebasic);
colorbar;
subplot(1, 2, 2);
imagesc(rerhr);
colorbar;