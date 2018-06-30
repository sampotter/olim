clear;

path(path, '../build/Release');

marks = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h'};

methodnames = {'basic', 'olim6_mp0', 'olim6_rhr', 'olim18_mp0', ...
               'olim18_rhr', 'olim26_mp0', 'olim26_rhr'};

ns = 2.^(2:6) + 1;
nmethods = length(methodnames);
T = zeros(length(ns), nmethods);

ntrials = 10;
for k = 1:length(ns)
    n = ns(k);
    fprintf('n = %d\n', n);
    for m = 1:nmethods
        methodname = methodnames{m};
        fprintf('- %s, trial:', methodname);
        B = zeros(n, n, n, 'logical');
        i0 = (n + 1)/2;
        B(i0, i0, i0) = 1;
        h = 2/(n - 1);
        T(k, m) = inf;
        for trial = 1:ntrials
            fprintf(' %d', trial);
            tic;
            fmm(B, 'Method', methodname, 'h', h, 'x0', 1, 'y0', 1, ...
                'z0', 1);
            T(k, m) = min(toc, T(k, m));
        end
        fprintf('\n');
    end
end

getplotsymb = @(index) strcat('-', marks{index});

figure;
loglog(ns, T(:, 1), getplotsymb(1));
hold on;
for m = 2:nmethods
    loglog(ns, T(:, m), getplotsymb(m));
end
ylabel('Time (s.)');
xlabel('n (s.t. # nodes = n^3)');

legendnames = {};
for m = 1:nmethods
    legendnames{m} = strrep(methodnames{m}, '_', '\_');
end
legend(legendnames);