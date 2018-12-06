function [Q, R] = qr3x2(A)
    assert(all(size(A) == [3 2]));
    Q = zeros(3, 2);
    R = zeros(2, 2);
    R(1, 1) = norm(A(:, 1));
    Q(:, 1) = A(:, 1)/R(1, 1);
    R(1, 2) = Q(:, 1)'*A(:, 2);
    z = A(:, 2) - R(1, 2)*Q(:, 1);
    R(2, 2) = norm(z);
    Q(:, 2) = z/R(2, 2);
end
