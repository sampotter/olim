#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--build_type', type=str, default='Release')
    args = parser.parse_args()

################################################################################
# preliminaries

if __name__ == '__main__':
    BUILD_TYPE = args.build_type
else:
    BUILD_TYPE = 'Release'

import sys
sys.path.insert(0, '../build/%s' % BUILD_TYPE)
sys.path.insert(0, '../misc/py')

import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import pyolim as olim
import time

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

plt.ion()
plt.style.use('bmh')

def tic():
    global t0
    t0 = time.perf_counter()

def toc():
    global t0
    return time.perf_counter() - t0

################################################################################
# gather timings

N = np.concatenate([
    np.logspace(3, 6, 12, base=2, dtype=int, endpoint=False),
    np.logspace(6, 9, 10, base=2, dtype=int)])

Tb, To = [], []

for i, n in enumerate(N):

    h = 2/(n-1)
    i0 = n//2
    S = np.ones((n, n, n))

    tb, to = np.inf, np.inf

    ntrials = 10 if n < 120 else 3

    for _ in range(ntrials):

        o = olim.BasicMarcher3D(S, h)
        o.add_boundary_node(i0, i0, i0)

        tic()
        o.run()
        tb = min(tb, toc())

        o = olim.Olim6Rect(S, h)
        o.add_boundary_node(i0, i0, i0)

        tic()
        o.run()
        to = min(to, toc())

    Tb.append(tb)
    To.append(to)

    print('n = %d, tb = %g, to = %g, to/tb = %g' % (n, tb, to, to/tb))

Tb, To = np.array(Tb), np.array(To)

np.savez('speed_comparison.npz', N=N, Tb=Tb, To=To)
