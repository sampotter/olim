function W = plot_winding_number(F, ndiv)
    k = 2^ndiv + 1;
    l = round(sqrt(2)*k);
    xy0 = [linspace(0, 1, k); zeros(1, k)];
    xy1 = [linspace(1, 0, l); linspace(0, 1, l)];
    xy2 = [zeros(1, k); linspace(1, 0, k)];
    xy = [xy0(:, 1:(size(xy0, 2) - 1)) xy1(:, 1:(size(xy1, 2) - 1)) ...
          xy2];
    theta = @(xy) angle(dot(F(xy), [1 1i]));

    th = zeros(1, size(xy, 2));
    for i = 1:size(xy, 2)
        th(i) = theta(xy(:, i));
    end
    th = unwrap(th);
    % h = [ones(1, k - 1)/(k - 1) ones(1, l - 1)/(l - 1) ones(1, k)/(k - 1)];
    
    figure;
    subplot(1, 2, 1);
    plot(th, 'k', 'LineWidth', 2);
    xlabel('k');
    ylabel('\theta');
    xlim([1 length(th)]);
    subplot(1, 2, 2);
    m = length(th);
    dth = th(2:m) - th(1:m-1);
    plot(dth, 'k', 'LineWidth', 2);
    xlabel('k');
    ylabel('d\theta');
    xlim([1 m]);
    
    W = sum(dth)/(2*pi);
end
