import sys
sys.path.insert(0, '../build/Release')

import eikonal as eik
import itertools
import matplotlib.pyplot as plt
import numpy as np
import speedfuncs
import speedfuncs3d
import time

marchers2d = [
    eik.BasicMarcher,
    eik.Olim4Mid0,
    eik.Olim4Rect,
    eik.Olim8Mid0,
    eik.Olim8Mid1,
    eik.Olim8Rect
];

marchers3d = [
    eik.BasicMarcher3D,
    eik.Olim6Mid0,
    eik.Olim6Mid1,
    eik.Olim6Rect,
    eik.Olim18Mid0,
    eik.Olim18Mid1,
    eik.Olim18Rect,
    eik.Olim26Mid0,
    eik.Olim26Mid1,
    eik.Olim26Rect
]

for marcher in marchers2d:
    marcher.dim = 2

for marcher in marchers3d:
    marcher.dim = 3

olim6_marchers = [
    eik.Olim6Mid0,
    eik.Olim6Mid1,
    eik.Olim6Rect
]

olim18_marchers = [
    eik.Olim18Mid0,
    eik.Olim18Mid1,
    eik.Olim18Rect
]

olim26_marchers = [
    eik.Olim26Mid0,
    eik.Olim26Mid1,
    eik.Olim26Rect
]

mid0_marchers = [
    eik.Olim6Mid0,
    eik.Olim18Mid0,
    eik.Olim26Mid0
]

mid1_marchers = [
    eik.Olim6Mid1,
    eik.Olim18Mid1,
    eik.Olim26Mid1
]

rect_marchers = [
    eik.Olim6Rect,
    eik.Olim18Rect,
    eik.Olim26Rect
]

_marcher_names = {
    eik.BasicMarcher3D: 'basic 3d',
    eik.BasicMarcher: 'basic marcher',
    eik.Olim18Mid0: 'olim18 mp0',
    eik.Olim18Mid1: 'olim18 mp1',
    eik.Olim18Rect: 'olim18 rhr',
    eik.Olim26Mid0: 'olim26 mp0',
    eik.Olim26Mid1: 'olim26 mp1',
    eik.Olim26Rect: 'olim26 rhr',
    eik.Olim4Mid0: 'olim4 mid0',
    eik.Olim4Rect: 'olim4 rect',
    eik.Olim6Mid0: 'olim6 mp0',
    eik.Olim6Mid1: 'olim6 mp1',
    eik.Olim6Rect: 'olim6 rhr',
    eik.Olim8Mid0: 'olim8 mid0',
    eik.Olim8Mid1: 'olim8 mid1',
    eik.Olim8Rect: 'olim8 rect',
}

def get_marcher_name(marcher):
    return _marcher_names[marcher]

def relerr(x, y, ord_):
    norm = lambda x: np.linalg.norm(x.flat, ord_)
    distxy = norm(x - y)
    return max(distxy/norm(x), distxy/norm(y))

def get_exact_soln(f, n):
    assert(f.dim in [2, 3])
    lin = np.linspace(-1, 1, n)
    if f.dim == 2:
        x, y = np.meshgrid(lin, lin)
        return f(x, y)
    else:
        x, y, z = np.meshgrid(lin, lin, lin)
        return f(x, y, z)

def compute_soln(marcher, s, n):
    i = int(n/2)
    h = 2/(n - 1)
    if marcher.dim == 2:
        m = marcher(n, n, h, s=s, x0=1, y0=1)
        m.addBoundaryNode(i, i)
    else:
        m = marcher(n, n, n, h, s=s, x0=1, y0=1, z0=1)
        m.addBoundaryNode(i, i, i)
    m.run()
    return np.array(m)

_line_style_cycler = itertools.cycle(['o-', 'o-.', 'o--', 'o:'])
def get_next_line_style():
    return next(_line_style_cycler);

def make_error_plot(s=speedfuncs.s1, minpow=3, maxpow=10,
                    marchers=marchers2d, verbose=True):
    assert(all(s.dim == marcher.dim for marcher in marchers))
    if verbose:
        print('make_error_plot:')
    f = speedfuncs.get_soln_func(s)
    ns = np.power(2, np.arange(minpow, maxpow + 1)) + 1
    E = {marcher: np.zeros(ns.shape) for marcher in marchers}
    for marcher, (i, n) in itertools.product(marchers, enumerate(ns)):
        u = get_exact_soln(f, n)
        U = compute_soln(marcher, s, n)
        e = relerr(u, U, np.inf)
        E[marcher][i] = e
        if verbose:
            print('- %s (n = %d, e = %g)' % (get_marcher_name(marcher), n, e))
    fig = plt.figure()
    for marcher in marchers:
        style = get_next_line_style()
        plt.loglog(ns, E[marcher], style, label=get_marcher_name(marcher))
    plt.legend()
    return fig

def tic():
    tic.t0 = time.time()
tic.t0 = None

def toc():
    if tic.t0:
        return time.time() - tic.t0
    else:
        raise RuntimeError("tic() hasn't been called")

def time_marcher(s, marcher, n, ntrials):
    assert(marcher.dim in [2, 3])
    t = np.inf
    h = 2/(n - 1)
    if marcher.dim == 2:
        for trial in range(ntrials):
            tic()
            marcher(n, n, h, x0=1, y0=1)
            t = min(t, toc())
    else:
        for trial in range(ntrials):
            tic()
            marcher(n, n, n, h, x0=1, y0=1, z0=1)
            t = min(t, toc())
    return t

def make_timing_plot(s=speedfuncs.s1, minpow=3, maxpow=10,
                     marchers=marchers2d, verbose=True, ntrials=10):
    assert(all(s.dim == marcher.dim for marcher in marchers))
    if verbose:
        print('make_timing_plot:')
    f = speedfuncs.get_soln_func(s)
    ns = np.power(2, np.arange(minpow, maxpow + 1)) + 1
    T = {marcher: np.zeros(ns.shape) for marcher in marchers}
    for marcher, (i, n) in itertools.product(marchers, enumerate(ns)):
        if verbose:
            print('- %s (n = %d)' % (get_marcher_name(marcher), n))
        T[marcher][i] = time_marcher(s, marcher, n, ntrials)
    fig = plt.figure()
    for marcher in marchers:
        style = get_next_line_style()
        plt.loglog(ns, T[marcher], style, label=get_marcher_name(marcher))
    plt.legend(loc=2)
    return fig

if __name__ == '__main__':
    # fig = make_error_plot(minpow=2, maxpow=7)
    fig = make_timing_plot(minpow=2, maxpow=12)
    fig.show()
