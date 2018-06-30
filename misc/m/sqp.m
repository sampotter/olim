function out = sqp(f, df, d2f, A, b, x0, tol, niters)
    n = size(A, 2);
    
    out = struct;
    out.xs = zeros(n, niters);
    out.xs(:, 1) = x0;
    out.dxs = zeros(n, niters);
    out.fs = zeros(1, niters);
    out.fs(1) = f(out.xs(:, 1));
    out.dfs = zeros(n, niters);
    out.alphas = zeros(1, niters);
    
    for k = 1:niters
        fprintf('k = %d\n', k);

        % Compute Hessian
        x = out.xs(:, k);
        H = d2f(x);

        % Perturb Hessian if it isn't positive definite
        lambda_min = min(eig(H));
        if lambda_min < 0
            fprintf('- fixing H\n');
            H = H - 1.1*lambda_min*eye(size(H));
            assert(min(eig(H)) > 0);
        end

        % Make Hessian symmetric (may only be necessary for quadprog)
        H = (H + H')/2;

        % Compute load vector for quadratic programx
        c = df(x) - H*x;

        % xopt = quadprog(H, c, A, b, [], [], [], [], zeros(N - 1, 1), opts);

        % Compute descent direction g
        found_opt = false;
        tol_ = eps;
        while ~found_opt
            try
                xopt = qpi(H, c, A, b, x, tol_, 10);
                found_opt = true;
            catch
                tol_ = 10*tol_;
                fprintf('- tol_ <- %g\n', tol_);
            end
            if tol_ > 1e-4
                error('some sort of failure happened')
            end
        end
        g = xopt - x;
        
        % Compute step size alpha
        c1 = 1e-4;
        alpha = 1;
        if norm(g, 'inf') > tol
            while f(x + alpha*g) > f(x) + c1*alpha*dot(df(x), g)
                alpha = 0.5*alpha;
                fprintf('- alpha = %g\n', alpha);
            end
        end

        out.alphas(k) = alpha;

        out.xs(:, k + 1) = x + alpha*g;
        out.dxs(:, k + 1) = out.xs(:, k + 1) - x;
        out.fs(k + 1) = f(out.xs(:, k + 1));
        out.dfs(:, k + 1) = out.fs(k + 1) - out.fs(k);
        
        if norm(out.dxs(:, k + 1), 'inf') < tol || ...
                abs(out.fs(k + 1) - out.fs(k)) < tol
            break
        end
    end

    out.iters = k;
end
