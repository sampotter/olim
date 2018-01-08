import sys

sys.path.insert(0, '../build/Release')

import eikonal as eik
import itertools
import matplotlib.pyplot as plt
import numpy as np
import speedfuncs

def relerr(x, y, ord_):
    norm = lambda x: np.linalg.norm(x, ord_)
    distxy = norm(x - y)
    return max(distxy/norm(x), distxy/norm(y))

def get_exact_soln(f, M):
    lin = np.linspace(-1, 1, M)
    x, y = np.meshgrid(lin, lin)
    return f(x, y)

def compute_soln(marcher, s, M):
    m = marcher(M, M, 2/(M - 1), s=s, x0=1, y0=1)
    i = int((M - 2)/2)
    m.addBoundaryNode(i, i)
    m.run()
    return np.array(m)

if __name__ == '__main__':
    s = speedfuncs.s1
    f = speedfuncs.f1

    minpow = 3
    maxpow = 10
    Ms = np.power(2, np.arange(minpow, maxpow + 1)) + 1

    marchers = [
        eik.BasicMarcher,
        eik.Olim4Mid0,
        eik.Olim4Rect,
        eik.Olim8Mid0,
        eik.Olim8Mid1,
        eik.Olim8Rect
    ]

    marcher_names = {
        eik.BasicMarcher: "basic marcher",
        eik.Olim4Mid0: "olim4 mid0",
        eik.Olim4Rect: "olim4 rect",
        eik.Olim8Mid0: "olim8 mid0",
        eik.Olim8Mid1: "olim8 mid1",
        eik.Olim8Rect: "olim8 rect"
    }

    E = {marcher: np.zeros(Ms.shape) for marcher in marchers}

    for marcher, (i, M) in itertools.product(marchers, enumerate(Ms)):
        u = get_exact_soln(f, M)
        U = compute_soln(marcher, s, M)
        e = relerr(u, U, np.inf)
        E[marcher][i] = e
        print("%s (M = %d, e = %g)" % (marcher_names[marcher], M, e))

    plt.figure()
    for marcher in marchers:
        plt.loglog(Ms, E[marcher], label=marcher_names[marcher])
    plt.legend()
    plt.show()
