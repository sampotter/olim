#!/usr/bin/env python3

BUILD_TYPE='Release'

import sys
if '../../build/%s' % BUILD_TYPE not in sys.path:
    sys.path.insert(0, '../../build/%s' % BUILD_TYPE)

import common
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs
import time

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

N = 2**np.arange(3, 9) + 1

use_local_factoring = True
r_fac = 0.01

Slows = [speedfuncs.s1, speedfuncs.s2, speedfuncs.s3, speedfuncs.s4]
Solns = {
    speedfuncs.s1: speedfuncs.f1,
    speedfuncs.s2: speedfuncs.f2,
    speedfuncs.s3: speedfuncs.f3,
    speedfuncs.s4: speedfuncs.f4
}
Olims = [eik.Olim4Mid0, eik.Olim4Mid1, eik.Olim4Rect,
         eik.Olim8Mid0, eik.Olim8Mid1, eik.Olim8Rect]

Slows_by_Olims = list(itertools.product(Slows, Olims))

T = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
E2 = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
EI = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}

ntrials = 3

for (slow, Olim), (ind, n) in itertools.product(Slows_by_Olims, enumerate(N)):

    print(Olim)

    # get timings

    h = 2/(n-1)
    i0, j0 = n//2, n//2
    L = np.linspace(-1, 1, n)
    x, y = np.meshgrid(L, L)
    R = np.sqrt(x**2 + y**2)
    I, J = np.where(R < r_fac)
    u, S = speedfuncs.get_fields(Solns[slow], slow, x, y)

    t = np.inf

    for _ in range(ntrials):

        o = Olim(S, h)
        if use_local_factoring:
            for i, j in zip(I, J):
                o.set_node_parent(i, j, i0, j0)
        o.addBoundaryNode(i0, j0)

        t0 = time.perf_counter()
        o.run()
        t = min(t, time.perf_counter() - t0)

    T[slow, Olim][ind] = t

    # get errors

    U = np.array([[o.getValue(i, j) for j in range(n)] for i in range(n)])
    E2[slow, Olim][ind] = norm(u - U, 'fro')/norm(u, 'fro')
    EI[slow, Olim][ind] = norm(u - U, np.inf)/norm(u, np.inf)

# make plots
 
fig, axes = plt.subplots(4, 2, sharex=True, figsize=(8, 8))

for row, slow in enumerate(Slows):
    for Olim in Olims:
        name = common.get_marcher_name(Olim)
        axes[row, 0].loglog(
            T[slow, Olim], E2[slow, Olim], '*-', linewidth=1, label=name)
        axes[row, 1].loglog(
            T[slow, Olim], EI[slow, Olim], '*-', linewidth=1, label=name)

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels)
fig.tight_layout()
fig.show()
