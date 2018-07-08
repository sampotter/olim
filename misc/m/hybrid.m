function [x, iters, deg] = hybrid(f, a, b, tol)
    if nargin < 4
        tol = eps;
    end
    
    deg = false;

    c = a;

    fa = f(a);
    fb = f(b);
    fc = f(c);
    
    assert(~(fa == 0 || fb == 0 || fc == 0));
    if sign(fb) == sign(fc)
        x = NaN;
        iters = [];
        deg = true;
        return;
    end

    iters = 0;
    while true
        if abs(fc) < abs(fb)
            tmp = c;  c = b;   b = tmp;  % swap(b, c)
            tmp = fc; fc = fb; fb = tmp; % swap(fb, fc)
            a = c;
            fa = fc;
        end
        iters = iters + 1;
        if abs(b - c) <= tol
            break;
        end
        dm = (c - b)/2;
        df = fa - fb;
        if df == 0
            ds = dm;
        else
            ds = -fb*(a - b)/df;
        end
        if sign(ds) ~= sign(dm) || abs(ds) > abs(dm)
            dd = dm;
        else
           dd = ds;
        end
        if abs(dd) < tol
            dd = tol*sign(dm)/2;
        end
        d = b + dd;
        fd = f(d);
        if fd == 0
            c = d;
            b = c;
            fc = fd;
            fb = fc;
            break;
        end
        a = b;
        b = d;
        fa = fb;
        fb = fd;
        if sign(fb) == sign(fc)
            c = a;
            fc = fa;
        end
    end
    x = (b + c)/2;
end
