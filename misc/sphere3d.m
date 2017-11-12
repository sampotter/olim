% A simple script that plots a sphere using each of the 3D fast
% marchers. This script is for use as a sanity check/debug script.

clear;

path(path, '../build/Release');

M = 5;
n = 2^M + 1;
B = zeros(n, n, n, 'logical');
B((n + 1)/2, (n + 1)/2, (n + 1)/2) = 1;
h = 2/(n - 1);
getU = @(method) fmm(B, 'Method', method, 'h', h, 'x0', 1, 'y0', 1, 'z0', 1);

[x y z] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, n), linspace(-1, 1, n));

% methodnames = {'basic', 'olim6_rhr', 'olim18_rhr', 'olim26_rhr'};
% methodnames = {'basic', 'olim6_mp0', 'olim18_mp0', 'olim26_mp0'};
% methodnames = {'basic', 'olim6_mp1', 'olim18_mp1', 'olim26_mp1'};

rows = 2;
cols = 2;
assert(length(methodnames) <= rows*cols);

for k = 1:length(methodnames)
    method = methodnames{k};
    U = getU(method);
    subplot(rows, cols, k);
    title(strrep(method, '_', '\_'));
    isosurface(x, y, z, U, 1);
    colormap gray;
    view([0 0]);
    camlight;
    lighting flat;
end
