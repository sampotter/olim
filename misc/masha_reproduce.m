path(path, '../build/Release');

B = zeros(21, 'logical');
B(11, 11) = 1;

h = 0.1;
s = @(x, y) masha_s(x, y);

Ubasic = fmm(B, 'h', h, 'Speed', s, 'Method', 'basic', 'x0', 1, 'y0', 1);

ubasic = importdata('gfun_basic.txt');
Ibasic = ubasic < 1e5; % mask region with valid solution
ubasic = circshift(ubasic, -1, 2); % align these two so that
Ibasic = circshift(Ibasic, -1, 2); % the minimum is at the center NB:
                                   % the minimum in the file is != 0,
                                   % but more like 1e-7... okay for
                                   % testing, though


tmp = zeros(21);
tmp(Ibasic) = abs(Ubasic(Ibasic) - ubasic(Ibasic));
imagesc(tmp)