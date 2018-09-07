#!/usr/bin/env python3

BUILD_TYPE='Release'

import os
import sys
if '../../build/%s' % BUILD_TYPE not in sys.path:
    sys.path.insert(0, os.path.abspath('../../build/%s' % BUILD_TYPE))

import common3d
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs3d
import time

from cycler import cycler
from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Npows = np.arange(3, 10)
N = 2**Npows + 1

use_local_factoring = True
r_fac = 0.1

Slows = [speedfuncs3d.s1, speedfuncs3d.s2, speedfuncs3d.s3, speedfuncs3d.s4]
Solns = {
    speedfuncs3d.s1: speedfuncs3d.f1,
    speedfuncs3d.s2: speedfuncs3d.f2,
    speedfuncs3d.s3: speedfuncs3d.f3,
    speedfuncs3d.s4: speedfuncs3d.f4
}
Olims = [eik.Olim6Mid0, eik.Olim6Mid1, eik.Olim6Rect,
         eik.Olim18Mid0, eik.Olim18Mid1, eik.Olim18Rect,
         eik.Olim26Mid0, eik.Olim26Mid1, eik.Olim26Rect,
         eik.Olim3dHuMid0, eik.Olim3dHuMid1, eik.Olim3dHuRect,
         eik.BasicMarcher3D]

Slows_by_Olims = list(itertools.product(Slows, Olims))

T = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
E2 = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
EI = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}

ntrials = 2

current_slow, current_Olim, current_n = None, None, None
for (slow, Olim), (ind, n) in itertools.product(Slows_by_Olims, enumerate(N)):
    if slow != current_slow:
        print(speedfuncs3d.get_slowness_func_name(slow))
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
    u, S = speedfuncs3d.get_fields(Solns[slow], slow, x, y, z)

    t = np.inf

    for _ in range(ntrials):

        o = Olim(S, h)
        if use_local_factoring:
            for i, j, k in zip(I, J, K):
                o.set_node_fac_parent(i, j, k, i0, j0, k0)
        o.addBoundaryNode(i0, j0, k0)

        t0 = time.perf_counter()
        o.run()
        t = min(t, time.perf_counter() - t0)

    T[slow, Olim][ind] = t

    # get errors

    U = np.array([[[o.getValue(i, j, k) for i in range(n)]
                   for j in range(n)]
                  for k in range(n)])
    E2[slow, Olim][ind] = norm((u - U).flatten())/norm(u.flatten())
    EI[slow, Olim][ind] = norm((u - U).flatten(), np.inf)/norm(u.flatten(), np.inf)

# make plots

marker = '*'
linestyles = ['solid', 'dashed', 'dotted']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

# time vs error

# make plots for the first two slowness functions (these are very
# cramped, so split them into two figures)

fig, axes = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(8, 8))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows[:2]):
    for ind, Olim in enumerate(Olims[:-1]):
        name = common3d.get_marcher_plot_name(Olim)
        axes[row, 0].loglog(
            T[slow, Olim], E2[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 0].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 1),
                          transform=axes[row, 0].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 1].loglog(
            T[slow, Olim], EI[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 1].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 1),
                          transform=axes[row, 1].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
    Olim = Olims[-1]
    axes[row, 0].loglog(T[slow, Olim], E2[slow, Olim], 'k-^', linewidth=1)
    axes[row, 1].loglog(T[slow, Olim], EI[slow, Olim], 'k-^', linewidth=1)

axes[-1, 0].set_xlabel('Time (s.)')    
axes[-1, 1].set_xlabel('Time (s.)')    

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=4)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.85)
fig.savefig('time_vs_error_3d_part1.eps')
fig.show()

# ... and the second two slowness functions

fig, axes = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(8, 8))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows[2:]):
    for ind, Olim in enumerate(Olims[:-1]):
        name = common3d.get_marcher_plot_name(Olim)
        axes[row, 0].loglog(
            T[slow, Olim], E2[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 0].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 3),
                          transform=axes[row, 0].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 1].loglog(
            T[slow, Olim], EI[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 1].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 3),
                          transform=axes[row, 1].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
    Olim = Olims[-1]
    axes[row, 0].loglog(T[slow, Olim], E2[slow, Olim], 'k-^', linewidth=1)
    axes[row, 1].loglog(T[slow, Olim], EI[slow, Olim], 'k-^', linewidth=1)

axes[-1, 0].set_xlabel('Time (s.)')    
axes[-1, 1].set_xlabel('Time (s.)')    

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=4)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.85)
fig.savefig('time_vs_error_3d_part2.eps')
fig.show()

# size vs error

# first two slowness functions...

fig, axes = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(8, 8))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows[:2]):
    for ind, Olim in enumerate(Olims[:-1]):
        name = common3d.get_marcher_plot_name(Olim)
        axes[row, 0].loglog(
            N, E2[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 0].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 1),
                          transform=axes[row, 0].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 0].minorticks_off()
        axes[row, 1].loglog(
            N, EI[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 1].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 1),
                          transform=axes[row, 1].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 1].minorticks_off()
    Olim = Olims[-1]
    axes[row, 0].loglog(N, E2[slow, Olim], 'k-^', linewidth=1)
    axes[row, 1].loglog(N, EI[slow, Olim], 'k-^', linewidth=1)

axes[-1, 0].set_xlabel('$N$')

xticklabels = ['$2^{%d} + 1$' % p for p in Npows]
axes[-1, 1].set_xlabel('$N$')    
axes[-1, 1].set_xticks(N)
axes[-1, 1].set_xticklabels(xticklabels)

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=4)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.85)
fig.savefig('size_vs_error_3d_part1.eps')
fig.show()

# ... and the second two slowness functions

fig, axes = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(8, 8))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows[2:]):
    for ind, Olim in enumerate(Olims[:-1]):
        name = common3d.get_marcher_plot_name(Olim)
        axes[row, 0].loglog(
            N, E2[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 0].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 3),
                          transform=axes[row, 0].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 0].minorticks_off()
        axes[row, 1].loglog(
            N, EI[slow, Olim], marker=marker, color=colors[ind//3],
            linestyle=linestyles[ind % 3], linewidth=1, label=name)
        axes[row, 1].text(0.95, 0.9, '$\\texttt{s%d}$' % (row + 3),
                          transform=axes[row, 1].transAxes,
                          horizontalalignment='center',
                          verticalalignment='center')
        axes[row, 1].minorticks_off()
    Olim = Olims[-1]
    axes[row, 0].loglog(N, E2[slow, Olim], 'k-^', linewidth=1)
    axes[row, 1].loglog(N, EI[slow, Olim], 'k-^', linewidth=1)

axes[-1, 0].set_xlabel('$N$')

xticklabels = ['$2^{%d} + 1$' % p for p in Npows]
axes[-1, 1].set_xlabel('$N$')    
axes[-1, 1].set_xticks(N)
axes[-1, 1].set_xticklabels(xticklabels)

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=4)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.85)
fig.savefig('size_vs_error_3d_part2.eps')
fig.show()
