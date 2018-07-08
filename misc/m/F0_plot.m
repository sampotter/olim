figure;

subplot(2, 2, 1);
if using_Z0_and_Z1
    hold on;
    contour(X, Y, Z0);
    contour(X, Y, Z1);
    contour(X, Y, Z0fac);
else
    contour(X, Y, Z); hold on;
end
plot([0 1 0 0], [0 0 1 0], 'k', 'LineWidth', 2); hold on;
if using_Z0_and_Z1
    plot(out0.xs(1, 1:out0.iters), out0.xs(2, 1:out0.iters), '-ok'); hold on;
    plot(out0.xs(1, out0.iters), out0.xs(2, out0.iters), '*k'); hold on;
    plot(out1.xs(1, 1:out1.iters), out1.xs(2, 1:out1.iters), '-om'); hold on;
    plot(out1.xs(1, out1.iters), out1.xs(2, out1.iters), '*m'); hold on;
    plot(out0fac.xs(1, 1:out0fac.iters), out0fac.xs(2, 1:out0fac.iters), '-og'); hold on;
    plot(out0fac.xs(1, out0fac.iters), out0fac.xs(2, out0fac.iters), '*g'); hold on;
else
    plot(out.xs(1, 1:out.iters), out.xs(2, 1:out.iters), '-o'); hold on;
    plot(out.xs(1, out.iters), out.xs(2, out.iters), '*m'); hold on;
end
colorbar;

subplot(2, 2, 2);
if using_Z0_and_Z1
    deltas0 = out0.fs(2:out0.iters) - out0.fs(1:(out0.iters - 1));
    semilogy(2:out0.iters, abs(deltas0), '*-k', 'LineWidth', 2);
    hold on;
    deltas1 = out1.fs(2:out1.iters) - out1.fs(1:(out1.iters - 1));
    semilogy(2:out1.iters, abs(deltas1), '*-m', 'LineWidth', 2);
    deltas0fac = out0fac.fs(2:out0fac.iters) - out0fac.fs(1:(out0fac.iters - 1));
    semilogy(2:out0fac.iters, abs(deltas0fac), '*-g', 'LineWidth', 2);
    legend('F0', 'F1', 'F0fac');
else
    deltas = out.fs(2:out.iters) - out.fs(1:(out.iters - 1));
    semilogy(2:out.iters, abs(deltas), '*-k', 'LineWidth', 2);
end
ylabel('dF_k');
xlabel('k');
if using_Z0_and_Z1
    xlim([1 max(out0.iters, out1.iters)]);
else
    xlim([1 out.iters]);
end

subplot(2, 2, 3);
if using_Z0_and_Z1
    plot(1:out0.iters, out0.fs(1:out0.iters), '*-k', 'LineWidth', 2);
    hold on;
    plot(1:out1.iters, out1.fs(1:out1.iters), '*-m', 'LineWidth', 2);
    plot(1:out0fac.iters, out0fac.fs(1:out0fac.iters), '*-g', 'LineWidth', 2);
    legend('F0', 'F1', 'F0fac');
else
    plot(1:out.iters, out.fs(1:out.iters), '*-k', 'LineWidth', 2);
end
ylabel('F_k');
xlabel('k');
if using_Z0_and_Z1
    xlim([1 max(out0.iters, out1.iters)]);
else
    xlim([1 out.iters]);
end

subplot(2, 2, 4);
if using_Z0_and_Z1
    plot(1:out0.iters, out0.alphas(1:out0.iters), '*-k', 'LineWidth', 2);
    hold on;
    plot(1:out1.iters, out1.alphas(1:out1.iters), '*-m', 'LineWidth', 2);
    plot(1:out0fac.iters, out0fac.alphas(1:out0fac.iters), '*-g', 'LineWidth', 2);
    legend('F0', 'F1', 'F0fac');
else
    plot(1:out.iters, out.alphas(1:out.iters), '*-k', 'LineWidth', 2);
end
ylabel('\alpha_k')
xlabel('k');
ylim([0.0 1.1]);
if using_Z0_and_Z1
    xlim([1 max(out0.iters, out1.iters)]);
else
    xlim([1 out.iters]);
end
