clear;

path(path, '../build/Release');

ntrials = 10;
mnum = 1;
Ms = 2*ceil(logspace(1, 3, 10)/2) + 1;
methods = {'basic', 'olim8_rhr', 'olim8_mp0', 'olim8_mp1'};
T = inf(length(Ms), length(methods));
for method = methods
    method = method{1};
    fprintf('method = %s\n', method);
    k = 1;
    for M = Ms
        fprintf('- M = %d\n', M);
        B = zeros(M, 'logical');
        B(1) = 1;
        for trial = 1:ntrials
            tic;
            fmm(B, 'Method', method);
            T(k, mnum) = min(T(k, mnum), toc);
        end
        k = k + 1;
    end
    mnum = mnum + 1;
end

figure;
loglog(Ms, T(:, 1)); hold on;
loglog(Ms, T(:, 2));
loglog(Ms, T(:, 3));
loglog(Ms, T(:, 4));
legend('basic', 'olim8\_rhr', 'olim8\_mp0', 'olim8\_mp1');