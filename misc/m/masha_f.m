function F = masha_f(X, Y)
    X = X - 0.900367222589747;
    F = (Y.^2 + (X + cos((X + Y)/2)).^2)/2;
end
