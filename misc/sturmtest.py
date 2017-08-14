import numpy as np

a = np.array([0.0864586, 0.0376806, -0.0385067, -0.00665833, 0.0036537])
b = np.polyder(a)
c = -np.polydiv(a, b)[1]
d = -np.polydiv(b, c)[1]
e = -np.polydiv(c, d)[1]

def sigma(x, polys):
    sign = np.sign(np.polyval(polys[0], x))
    changes = 0
    for i in range(1, len(polys)):
        newsign = np.sign(np.polyval(polys[i], x))
        if sign != 0 and newsign != 0 and sign != newsign:
            changes += 1
        sign = newsign
    return changes

def sturm(a, b, polys):
    assert(np.polyval(polys[0], a) != 0)
    assert(np.polyval(polys[0], b) != 0)
    return sigma(a, polys) - sigma(b, polys)

def secant(f, x0, x1, a=0, b=1, tol=1e-13):
    x = (x1*f(x0) - x0*f(x1))/(f(x0) - f(x1))
    x1 = x0
    x0 = x
    if x < a or b < x:
        return x, abs(f(x)) <= tol
    while abs(f(x)) > tol:
        x = (x1*f(x0) - x0*f(x1))/(f(x0) - f(x1))
        x1 = x0
        x0 = x
        if x < a or b < x:
            return x, abs(f(x)) <= tol
    return x, True

def dst(poly):
    pass

def findroots(poly, a, b):
    polys = [poly]
    polys.append(np.polyder(poly))
    for i in range(3):
        polys.append(-np.polydiv(polys[i], polys[i + 1])[1])

    def f(x): return np.polyval(poly, x)
        
    def rec(a, b):
        nroots = sturm(a, b, polys)
        if nroots == 0:
            return []
        elif nroots == 1:
            h = 0.1
            foundroot = False
            root, foundroot = secant(f, a, a + (b - a)*h, a, b)
            if not foundroot:
                root, foundroot = secant(f, b, b + (a - b)*h, a, b)
            if not foundroot:
                raise Exception("couldn't find root: a = %g, b = %g, h = %g" % (a, b, h))
            return [root]
        else:
            c = (a + b)/2
            return rec(a, c) + rec(c, b)

    return rec(a, b)
