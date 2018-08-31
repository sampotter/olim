import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import pyeikonal as eik
import numpy as np
import time

marchers = [eik.BasicMarcher, eik.Olim4Mid0, eik.Olim4Mid1, eik.Olim4Rect,
            eik.Olim8Mid0, eik.Olim8Mid1, eik.Olim8Rect]

olim4_marchers = [eik.Olim4Mid0, eik.Olim4Mid1, eik.Olim4Rect]
olim8_marchers = [eik.Olim8Mid0, eik.Olim8Mid1, eik.Olim8Rect]

mid0_marchers = [eik.Olim4Mid0, eik.Olim8Mid0]
mid1_marchers = [eik.Olim4Mid1, eik.Olim8Mid1]
rect_marchers = [eik.Olim4Rect, eik.Olim8Rect]

_marcher_names = {
    eik.BasicMarcher: 'basic',
    eik.Olim4Mid0: 'olim4 mp0',
    eik.Olim4Mid1: 'olim4 mp1',
    eik.Olim4Rect: 'olim4 rhr',
    eik.Olim8Mid0: 'olim8 mp0',
    eik.Olim8Mid1: 'olim8 mp1',
    eik.Olim8Rect: 'olim8 rhr'}

def get_marcher_name(marcher):
    return _marcher_names[marcher]

_marcher_plot_names = {
    eik.BasicMarcher: '\\texttt{basic}',
    eik.Olim4Mid0: '\\texttt{olim4\_mp0}',
    eik.Olim4Mid1: '\\texttt{olim4\_mp1}',
    eik.Olim4Rect: '\\texttt{olim4\_rhr}',
    eik.Olim8Mid0: '\\texttt{olim8\_mp0}',
    eik.Olim8Mid1: '\\texttt{olim8\_mp1}',
    eik.Olim8Rect: '\\texttt{olim8\_rhr}'}

def get_marcher_plot_name(marcher):
    return _marcher_plot_names[marcher]

_marchers_by_name = {v: k for k, v in _marcher_names.items()}

def get_marcher_by_name(name):
    return _marchers_by_name[name]

def relerr(x, y, ord_):
    norm = lambda x: np.linalg.norm(x.flat, ord_)
    distxy = norm(x - y)
    return max(distxy/norm(x), distxy/norm(y))

def get_exact_soln(f, M):
    l = np.linspace(-1, 1, M)
    return f(*np.meshgrid(l, l))

def compute_soln(marcher, s, M):
    l = np.linspace(-1, 1, M)
    m = marcher(s(*np.meshgrid(l, l)), 2/(M - 1))
    m.addBoundaryNode(int(M/2), int(M/2))
    m.run()
    U = np.array([[m.getValue(i, j) for j in range(M)] for i in range(M)])
    return U

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
    x, y = np.meshgrid(l, l)
    s_cache = s(x, y)
    def do_trial(trial):
        tic()
        m = Marcher(s_cache, h)
        m.addBoundaryNode(i, i)
        m.run()
        return toc()
    times = min(do_trial(t) for t in range(ntrials))
    return times
