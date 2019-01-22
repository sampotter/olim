import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import pyolim as olim
import numpy as np
import time

marchers = [
    olim.BasicMarcher3D,
    olim.Olim6Mid0,
    olim.Olim6Mid1,
    olim.Olim6Rect,
    olim.Olim18Mid0,
    olim.Olim18Mid1,
    olim.Olim18Rect,
    olim.Olim26Mid0,
    olim.Olim26Mid1,
    olim.Olim26Rect,
    olim.Olim3dHuMid0,
    olim.Olim3dHuMid1,
    olim.Olim3dHuRect
]

olim6_marchers = [
    olim.Olim6Mid0,
    olim.Olim6Mid1,
    olim.Olim6Rect
]

olim18_marchers = [
    olim.Olim18Mid0,
    olim.Olim18Mid1,
    olim.Olim18Rect
]

olim26_marchers = [
    olim.Olim26Mid0,
    olim.Olim26Mid1,
    olim.Olim26Rect
]

mid0_marchers = [
    olim.Olim6Mid0,
    olim.Olim18Mid0,
    olim.Olim26Mid0
]

mid1_marchers = [
    olim.Olim6Mid1,
    olim.Olim18Mid1,
    olim.Olim26Mid1
]

rect_marchers = [
    olim.Olim6Rect,
    olim.Olim18Rect,
    olim.Olim26Rect
]

_marcher_names = {
    olim.BasicMarcher3D: 'basic 3d',
    olim.Olim6Mid0: 'olim6 mp0',
    olim.Olim6Mid1: 'olim6 mp1',
    olim.Olim6Rect: 'olim6 rhr',
    olim.Olim18Mid0: 'olim18 mp0',
    olim.Olim18Mid1: 'olim18 mp1',
    olim.Olim18Rect: 'olim18 rhr',
    olim.Olim26Mid0: 'olim26 mp0',
    olim.Olim26Mid1: 'olim26 mp1',
    olim.Olim26Rect: 'olim26 rhr',
    olim.Olim3dHuMid0: 'olim3d hu mp0',
    olim.Olim3dHuMid1: 'olim3d hu mp1',
    olim.Olim3dHuRect: 'olim3d hu rhr'
}

def get_marcher_name(marcher):
    return _marcher_names[marcher]

_marchers_by_name = {v: k for k, v in _marcher_names.items()}

def get_marcher_by_name(name):
    return _marchers_by_name[name]

_marcher_plot_names = {
    olim.BasicMarcher3D: '\\texttt{basic\_3d}',
    olim.Olim6Mid0: '\\texttt{olim6\_mp0}',
    olim.Olim6Mid1: '\\texttt{olim6\_mp1}',
    olim.Olim6Rect: '\\texttt{olim6\_rhr}',
    olim.Olim18Mid0: '\\texttt{olim18\_mp0}',
    olim.Olim18Mid1: '\\texttt{olim18\_mp1}',
    olim.Olim18Rect: '\\texttt{olim18\_rhr}',
    olim.Olim26Mid0: '\\texttt{olim26\_mp0}',
    olim.Olim26Mid1: '\\texttt{olim26\_mp1}',
    olim.Olim26Rect: '\\texttt{olim26\_rhr}',
    olim.Olim3dHuMid0: '\\texttt{olim3d\_mp0}',
    olim.Olim3dHuMid1: '\\texttt{olim3d\_mp1}',
    olim.Olim3dHuRect: '\\texttt{olim3d\_rhr}'
}

def get_marcher_plot_name(marcher):
    return _marcher_plot_names[marcher]

def relerr(x, y, ord_):
    norm = lambda x: np.linalg.norm(x.flat, ord_)
    distxy = norm(x - y)
    return max(distxy/norm(x), distxy/norm(y))

def get_exact_soln(f, M):
    l = np.linspace(-1, 1, M)
    return f(*np.meshgrid(l, l, l))

def compute_soln(marcher, s, M):
    l = np.linspace(-1, 1, M)
    m = marcher(s(*np.meshgrid(l, l, l)), 2/(M - 1))
    m.addBoundaryNode(int(M/2), int(M/2), int(M/2))
    m.run()
    return np.array(m)

def tic():
    tic.t0 = time.time()
tic.t0 = None

def toc():
    if tic.t0:
        return time.time() - tic.t0
    else:
        raise RuntimeError("tic() hasn't been called")

def time_marcher(Marcher, s, n, ntrials=10):
    print('  - n = %d' % n)
    h = 2/(n - 1)
    l = np.linspace(-1, 1, n)
    i = int(n/2)
    x, y, z = np.meshgrid(l, l, l)
    s_cache = s(x, y, z)
    def do_trial(trial):
        tic()
        m = Marcher(s_cache, h)
        m.addBoundaryNode(i, i, i)
        m.run()
        return toc()
    times = min(do_trial(t) for t in range(ntrials))
    return times
