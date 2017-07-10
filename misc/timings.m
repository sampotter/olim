clear;

path(path, '../build/Release');

s1 = @(x, y) 1 - sin(sqrt(x.^2 + y.^2));

ntrials = 10;
mnum = 1;
Ms = 2*ceil(logspace(1, 3, 10)/2) + 1;
methods = {'basic', 'olim4_rhr', 'olim4_rhr_lut', 'olim4_mp0', ...
           'olim8_rhr', 'olim8_mp0', 'olim8_mp1'};
T = inf(length(Ms), length(methods));
for k = 1:length(methods)
    method = methods{k};
    fprintf('method = %s\n', method);
    k = 1;
    for M = Ms
        fprintf('- M = %d\n', M);
        B = zeros(M, 'logical');
        B((M + 1)/2, (M + 1)/2) = 1;
        h = 2/(M - 1);
        for trial = 1:ntrials
            tic;
            fmm(B, 'h', h, 'Method', method, 'Speed', s1, 'x0', 1, 'y0', 1);
            T(k, mnum) = min(T(k, mnum), toc);
        end
        k = k + 1;
    end
    mnum = mnum + 1;
end

figure;
for k = 1:length(methods)
    loglog(Ms, T(:, k));
    hold on;
end
legend('basic', 'olim4\_rhr', 'olim4\_lut\_rhr', 'olim4\_mp0', ...
       'olim8\_rhr', 'olim8\_mp0', 'olim8\_mp1');
