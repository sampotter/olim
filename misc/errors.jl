include("../src/fmm.jl")

s1(x, y) = 1 - sin(sqrt(x.^2 + y.^2));
f1(x, y) = cos(sqrt(x.^2 + y.^2)) + sqrt(x.^2 + y.^2) - 1;

s = s1;
f = f1;

Ms = (2.^(3:10)) + 1;
for M = Ms
    println("M = ", M)

    B = Array{Bool, 2}(zeros(M, M))
    B[Int((M + 1)/2), Int((M + 1)/2)] = true;

    h = 2/(M - 1);

    u = zeros(M, M);
    xs = linspace(-1, 1, M)
    ys = linspace(-1, 1, M)
    for (i, y) = enumerate(ys), (j, x) = enumerate(xs)
        u[i, j] = f(x, y)
    end

    U = fmm.march(B, h=h, speed=s, x0=1, y0=1)
end
