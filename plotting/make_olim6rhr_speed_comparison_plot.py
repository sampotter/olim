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
import pyeikonal as eik
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

Nmin, Nmax = 11, 401
N = np.arange(Nmin, Nmax, 2)

ntrials = 3

Tb, To = [], []

for i, n in enumerate(N):

    h = 2/(n-1)
    i0 = n//2
    S = np.ones((n, n, n))

    tb, to = np.inf, np.inf

    for _ in range(ntrials):

        o = eik.BasicMarcher3D(S, h)
        o.addBoundaryNode(i0, i0, i0)

        tic()
        o.run()
        tb = min(tb, toc())

        o = eik.Olim6Rect(S, h)
        o.addBoundaryNode(i0, i0, i0)

        tic()
        o.run()
        to = min(to, toc())

    print('n = %d, tb = %g, to = %g, to/tb = %g' % (n, tb, to, to/tb))

Tb, To = np.array(Tb), np.array(To)

# TODO: also plot To and Tb to show effect of cache spill?

np.savez('speed_comparison.npz', N=N, Tb=Tb, To=To)

# plt.figure()
# plt.semilogx(N, To/Tb)
# plt.savefig('olim6rhr_speed_comparison_plot.eps')


# N = [9, 11, 13, 15, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 73, 81,
#      89, 97, 105, 113, 121, 129, 145, 161, 177, 193, 209, 225, 241, 257, 321,
#      385]

# R = [1.3441, 1.58329, 1.42205, 2.54741, 2.12998, 1.80207, 1.59377, 1.43226,
#      1.76175, 1.55943, 1.5695, 1.60596, 1.48476, 1.30242, 1.4271, 1.34489,
#      1.41143, 1.2998, 1.29134, 1.30058, 1.30324, 1.27189, 1.2581, 1.23572,
#      1.27774, 1.23627, 1.20808, 1.19071, 1.10093, 0.742436, 1.15878, 1.13266,
#      1.12554, 1.12279, 1.11564]

# plt.figure()
# plt.semilogx(N, R)
# plt.show()
