#!/usr/bin/env python3


################################################################################
# parse arguments first

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--min_3d_power', type=int, default=3)
parser.add_argument('--max_3d_power', type=int, default=10)
parser.add_argument('--build_type', type=str, default='Release')
parser.add_argument('--no_factoring', action='store_true')
args = parser.parse_args()

################################################################################
# preliminaries

import sys
sys.path.insert(0, '../build/%s' % args.build_type)
sys.path.insert(0, '../misc/py')

import common3d
import datetime
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyolim as olim
import slow3
import time

from cycler import cycler
from matplotlib import rc

rc('text', usetex=True)
rc('font', **{
    'family': 'serif',
    'serif': ['Computer Modern'],
    'size': 8
})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Npows = np.arange(args.min_3d_power, args.max_3d_power + 1)
N = 2**Npows + 1

r_fac = 0.1

Slows = [slow3.s1, slow3.s2, slow3.s3, slow3.s4]
Solns = {
    slow3.s1: slow3.f1,
    slow3.s2: slow3.f2,
    slow3.s3: slow3.f3,
    slow3.s4: slow3.f4
}
Olims = [olim.Olim6Mid0, olim.Olim6Mid1, olim.Olim6Rect,
         olim.Olim18Mid0, olim.Olim18Mid1, olim.Olim18Rect,
         olim.Olim26Mid0, olim.Olim26Mid1, olim.Olim26Rect,
         olim.Olim3dHuMid0, olim.Olim3dHuMid1, olim.Olim3dHuRect]

Slows_by_Olims = list(itertools.product(Slows, Olims))

T = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
E = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}

ntrials = 2

current_slow, current_Olim, current_n = None, None, None
for (slow, Olim), (ind, n) in itertools.product(Slows_by_Olims, enumerate(N)):
    if slow != current_slow:
        print(slow3.get_slowness_func_name(slow))
        current_slow = slow
    if Olim != current_Olim:
        print('* %s' % str(Olim))
        current_Olim = Olim
    if n != current_n:
        print('  - %d' % n)
        current_n = n

    # get timings

    h = 2/(n-1)
    i0, j0, k0 = n//2, n//2, n//2
    L = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(L, L, L)
    R = np.sqrt(x**2 + y**2 + z**2)
    I, J, K = np.where(R < r_fac)
    u, S = slow3.get_fields(Solns[slow], slow, x, y, z)

    t = np.inf

    for _ in range(1 if n > 100 else ntrials):

        o = Olim(S, h)
        if not args.no_factoring:
            fc = olim.FacCenter3d(i0, j0, k0, slow(0, 0, 0))
            for i, j, k in zip(I, J, K):
                o.set_fac_src(i, j, k, fc)
        o.add_boundary_node(i0, j0, k0)

        t0 = time.perf_counter()
        o.run()
        t = min(t, time.perf_counter() - t0)

        print('    + %s' % datetime.timedelta(seconds=t))

    T[slow, Olim][ind] = t

    # get errors

    U = np.array([[[o.get_value(i, j, k) for k in range(n)]
                   for j in range(n)]
                  for i in range(n)])
    E[slow, Olim][ind] = \
        norm((u - U).flatten(), np.inf)/norm(u.flatten(), np.inf)

# make plots

marker = '|'
linestyles = ['solid', 'dashed', 'dotted']
# colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
cmap = [0, 1, 4, 3]
# linestyles = ['-', '--', ':', '-.']
linestyles = [':', '-.', '--', '-']

# time vs error

# make plots for the first two slowness functions (these are very
# cramped, so split them into two figures)

fig, axes = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(6.5, 6))

ax = axes.flatten()

for i, slow in enumerate(Slows):
    for ind, Olim in enumerate(Olims):
        print((i, ind, Olim))
        ax[i].loglog(
            T[slow, Olim], E[slow, Olim], marker=marker, markersize=3.5,
            color=colors[cmap[ind % 3]], linestyle=linestyles[ind//3],
            linewidth=1, label=common3d.get_marcher_plot_name(Olim))
        ax[i].text(
            0.95, 0.9, '$\\texttt{s%d}$' % (i + 1),
            transform=ax[i].transAxes, horizontalalignment='center',
            verticalalignment='center')

axes[0, 0].set_ylabel(r'$\|u - U\|_\infty/\|u\|_\infty$')
axes[1, 0].set_ylabel(r'$\|u - U\|_\infty/\|u\|_\infty$')
axes[-1, 0].set_xlabel('Time (s.)')    
axes[-1, 1].set_xlabel('Time (s.)')    

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=4)
fig.tight_layout()
fig.subplots_adjust(0.085, 0.055, 0.995, 0.935)
fig.show()

fig.savefig('time_vs_error_3d.eps')
