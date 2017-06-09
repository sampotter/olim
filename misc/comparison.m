B = zeros(201, 'logical'); 
B(101, 101) = 1;
h = 0.01;
F = @(x, y) 1;
U1 = fmm(B, h, F, 'basic');
U2 = fmm(B, h, F, 'olim8pt');
[x y] = meshgrid(linspace(-1, 1, 201), linspace(-1, 1, 201));
u = sqrt(x.^2 + y.^2);
rms(u(:) - U1(:))
rms(u(:) - U2(:))
norm(u(:) - U1(:), 'inf')
norm(u(:) - U2(:), 'inf')