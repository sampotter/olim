figure;

subplot(2, 2, 1);
contour(X, Y, Z); hold on;
plot([0 1 0 0], [0 0 1 0], 'k', 'LineWidth', 2); hold on;
plot(out.xs(1, 1:out.iters), out.xs(2, 1:out.iters), '-o'); hold on;
plot(out.xs(1, out.iters), out.xs(2, out.iters), '*m'); hold on;
colorbar;

subplot(2, 2, 2);
deltas = out.fs(2:out.iters) - out.fs(1:(out.iters - 1));
semilogy(2:out.iters, abs(deltas), '*-k', 'LineWidth', 2);
ylabel('dF_k');
xlabel('k');
xlim([1 out.iters]);

subplot(2, 2, 3);
plot(1:out.iters, out.fs(1:out.iters), '*-k', 'LineWidth', 2);
ylabel('F_k');
xlabel('k');
xlim([1 out.iters]);

subplot(2, 2, 4);
plot(1:out.iters, out.alphas(1:out.iters), '*-k', 'LineWidth', 2);
ylabel('\alpha_k')
xlabel('k');
ylim([0.0 1.1]);
xlim([1 out.iters]);