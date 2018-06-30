function make_s10_gt(M)
    path(path, '../build/Release'); % TODO: only do this if the
                                    % path doesn't already contain it

    B = zeros(M, 'logical');
    B((M + 1)/2, (M + 1)/2) = 1;
    h = 2/(M - 1);
    s = @(x, y) star_s(x, y);
    U = fmm(B, 'h', h, 'Speed', s, 'Method', 'olim8_rhr', 'x0', 1, 'y0', 1);
    save('U_s10_gt.mat', 'M', 'U');
end
