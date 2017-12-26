function x = qpi(H, c, A, b, x0, tol, niters)
% This is an implementation of the active set method for solving
% a quadratic program with linear inequality constraints (following
% Nocedal & Wright).
% 
% In this implementation, it's assumed that ALL of the constraints are
% inequality constraints, so the index set \mathcal{I} referred to in
% Nocedal & Wright is redundant and can be ignored.

    assert(issymmetric(H));
    n = size(H, 1);
    
    m = size(A, 1);

    if nargin <= 4
        x0 = zeros(n, 1);
    end
    x = x0;
    
    if nargin <= 5
        tol = eps;
    end
    
    if nargin <= 6
        niters = inf;
    end

    W = find(abs(A*x - b) <= tol);

    k = 1;
    while true
        p = qpe(H, H*x + c, A(W, :), zeros(length(W), 1));
        if norm(p, 'inf') <= tol
            lam = nan(m, 1);
            lam(W) = A(W, :)'\(H*x + c);
            if all(lam(W) >= 0)
                break
            end
            [~, j] = min(lam);
            W = setdiff(W, j);
        else
            V = setdiff(1:m, W);
            alpha = 1;
            denom = nan(m, 1);
            denom(V) = A(V, :)*p;
            ind = denom < 0;
            if any(ind)
                numer = nan(m, 1);
                numer(V) = b(V) - A(V, :)*x;
                mask = zeros(m, 1);
                mask(~ind) = nan;
                [minval, argmin] = min(mask + numer./denom);
                alpha = max(0, min(minval, alpha));
                if alpha < 1
                    W = union(W, argmin);
                end
            end
            x = x + alpha*p;
        end
        if k == niters
            error(sprintf('Failed to converge in %d steps', niters))
        end
        k = k + 1;
    end
end
