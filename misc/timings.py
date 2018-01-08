import sys

sys.path.insert(0, '../build/Release')

import eikonal as eik
import matplotlib.pyplot as plt
import numpy as np
import itertools
import speedfuncs
import time

def tic():
    tic.t0 = time.time()
tic.t0 = None

def toc():
    if tic.t0:
        return time.time() - tic.t0
    else:
        raise RuntimeError("tic() hasn't been called")

if __name__ == '__main__':
    s = speedfuncs.s2
    f = speedfuncs.f2

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

    ntrials = 10
    T = {marcher: np.zeros(Ms.shape) for marcher in marchers}

    for marcher, (i, M) in itertools.product(marchers, enumerate(Ms)):
        print('%s (M = %d)' % (marcher_names[marcher], M))
        t = np.inf
        for trial in range(ntrials):
            tic()
            marcher(M, M, 2/(M - 1), s, x0=1, y0=1)
            t = min(t, toc())
        T[marcher][i] = t

    plt.figure()
    for marcher in marchers:
        plt.loglog(Ms, T[marcher], label=marcher_names[marcher])
    plt.legend()
    plt.show()
