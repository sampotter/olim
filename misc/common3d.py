import eikonal as eik
import numpy as np
import time

marchers = [
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
    eik.Olim6Mid0: 'olim6 mp0',
    eik.Olim6Mid1: 'olim6 mp1',
    eik.Olim6Rect: 'olim6 rhr',
    eik.Olim18Mid0: 'olim18 mp0',
    eik.Olim18Mid1: 'olim18 mp1',
    eik.Olim18Rect: 'olim18 rhr',
    eik.Olim26Mid0: 'olim26 mp0',
    eik.Olim26Mid1: 'olim26 mp1',
    eik.Olim26Rect: 'olim26 rhr'
}

def get_marcher_name(marcher):
    return _marcher_names[marcher]

_marchers_by_name = {v: k for k, v in _marcher_names.items()}

def get_marcher_by_name(name):
    return _marchers_by_name[name]

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
    h = 2/(n - 1)
    l = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(l, l, l)
    s_cache = s(x, y, z)
    def do_trial():
        tic()
        Marcher(s_cache, h)
        return toc()
    return min(do_trial() for _ in range(ntrials))
