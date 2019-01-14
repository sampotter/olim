for i = 1:3
    for j = 1:3
        fprintf('p%d[%d] = %0.16g;\n', i-1, j-1, P(j, i));
    end
end
fprintf('u0 = %0.16g;\n', u0);
fprintf('u1 = %0.16g;\n', u1);
fprintf('u2 = %0.16g;\n', u2);
fprintf('s = %0.16g;\n', s);
fprintf('s0 = %0.16g;\n', s0);
fprintf('s1 = %0.16g;\n', s1);
fprintf('s2 = %0.16g;\n', s2);
fprintf('h = %0.16g;\n', h);
for i = 1:2
    fprintf('info.lambda[%d] = %0.16g;\n', i-1, out0.xs(i, out0.iters + 1));
end