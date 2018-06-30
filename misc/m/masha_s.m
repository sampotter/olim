function S = masha_s(X, Y)
    X = X - 0.900367222589747;
    aux0 = (X + Y)/2;
    aux1 = X + cos(aux0);
    aux2 = sin(aux0)/2;
    dX = aux1.*(1 - aux2);
    dY = Y - aux1.*aux2;
    S = sqrt(dX.^2 + dY.^2);
end
