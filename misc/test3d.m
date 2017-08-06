clear;

path(path, '../build/Release');

n = 31;
B = zeros(n, n, n, 'logical');
bi = (n + 1)/2;
B(bi, bi, bi) = 1;
L = -(n - 1)/2:(n - 1)/2;
[X Y Z] = meshgrid(L, L, L);
U = fmm(B);
isosurface(X, Y, Z, U, (n - 1)/2);
xlim([-bi, bi]);
ylim([-bi, bi]);
zlim([-bi, bi]);
lighting flat;
