function tetra222rhr(u0, u1, u2, s, h)
    path(path, 'conicsintersection/matlab');
    
    assert(nargin >= 3);
    if nargin <= 3, s = 1; end
    if nargin <= 4, h = 1; end

    a = abs(u1 - u0)/(s*h); a2 = a*a;
    b = abs(u2 - u0)/(s*h); b2 = b*b;

    Q1coefs = [2*a2 - 4; 2*a2 - 4; 2*a2 - 1; 4 - 2*a2; 2 - 2*a2; 2*a2 - 1];
    Q2coefs = [2*b2 - 1; 2*b2 - 4; 2*b2 - 4; 2 - 2*b2; 4 - 2*b2; 2*b2 - 1];

    ceval = @(x, y, C) C(1)*x.^2 + C(2)*x.*y + C(3)*y.^2 + C(4)*x + C(5)*y + C(6);
    cmatrix = @(C) [C(1) C(2)/2 C(4)/2; C(2)/2 C(3) C(5)/2; C(4)/2 C(5)/2 C(6)];

    P = intersectConics(cmatrix(Q1coefs), cmatrix(Q2coefs));
    P = [P(1, :)./P(3, :); P(2, :)./P(3, :)]
    Px = P(1, :)';
    Py = P(2, :)';

    L = linspace(-0.5, 1.5, 201);
    [x y] = meshgrid(L, L);

    figure;
    subplot(1, 1, 1);
    hold on;
    contour(x, y, ceval(x, y, Q1coefs));
    contour(x, y, ceval(x, y, Q2coefs));
    colorbar;
    plot([0, 1, 0, 0], [0, 0, 1, 0], 'k');
    scatter(Px, Py, 'k');
end
