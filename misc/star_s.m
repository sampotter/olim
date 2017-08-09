function s = star_s(x, y)
    x = 3*x;
    y = 3*y;
    z = x + i*y;
    r = abs(z);
    theta = angle(z);
    a = theta + sin(3*r).*cos(3*theta);
    b = r.*(1 + cos(3*a)/2);
    s = 1 + 2*exp(-b.*b/2);
end
