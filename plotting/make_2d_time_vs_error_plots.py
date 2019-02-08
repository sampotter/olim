#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--min_2d_power', type=int, default=3)
parser.add_argument('--max_2d_power', type=int, default=15)
parser.add_argument('--build_type', type=str, default='Release')
args = parser.parse_args()

################################################################################
# preliminaries

import sys
sys.path.insert(0, '../build/%s' % args.build_type)
sys.path.insert(0, '../misc/py')

import common
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyolim as olim
import slow2
import time

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

Npows = np.arange(args.min_2d_power, args.max_2d_power + 1)
N = 2**Npows + 1

use_local_factoring = True
r_fac = 0.1

Slows = [slow2.s1, slow2.s2, slow2.s3, slow2.s4]
Solns = {
    slow2.s1: slow2.f1,
    slow2.s2: slow2.f2,
    slow2.s3: slow2.f3,
    slow2.s4: slow2.f4
}
Olims = [olim.Olim4Mid0, olim.Olim4Mid1, olim.Olim4Rect,
         olim.Olim8Mid0, olim.Olim8Mid1, olim.Olim8Rect]

Slows_by_Olims = list(itertools.product(Slows, Olims))

T = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
E2 = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}
EI = {(slow, Olim): np.empty(N.shape) for slow, Olim in Slows_by_Olims}

ntrials = 2

current_slow, current_Olim, current_n = None, None, None
for (slow, Olim), (ind, n) in itertools.product(Slows_by_Olims, enumerate(N)):
    if slow != current_slow:
        print(slow2.get_slowness_func_name(slow))
        current_slow = slow
    if Olim != current_Olim:
        print('* %s' % str(Olim))
        current_Olim = Olim
    if n != current_n:
        print('  - %d' % n)
        current_n = n

    # get timings

    h = 2/(n-1)
    i0, j0 = n//2, n//2
    L = np.linspace(-1, 1, n)
    x, y = np.meshgrid(L, L)
    R = np.sqrt(x**2 + y**2)
    I, J = np.where(R < r_fac)
    u, S = slow2.get_fields(Solns[slow], slow, x, y)

    t = np.inf

    for _ in range(ntrials):

        o = Olim(S, h)
        if use_local_factoring:
            fc = olim.FacCenter(i0, j0, slow(0, 0))
            for i, j in zip(I, J):
                o.set_fac_src(i, j, fc)
        o.add_boundary_node(i0, j0)

        t0 = time.perf_counter()
        o.run()
        t = min(t, time.perf_counter() - t0)

    T[slow, Olim][ind] = t

    # get errors

    U = np.array([[o.get_value(i, j) for j in range(n)] for i in range(n)])
    E2[slow, Olim][ind] = norm(u - U, 'fro')/norm(u, 'fro')
    EI[slow, Olim][ind] = norm(u - U, np.inf)/norm(u, np.inf)

# make plots
 
marker = '*'
linestyles = ['solid', 'dashed', 'dotted']
colors = ['#1f77b4', '#d62728']

# time vs error

fig, axes = plt.subplots(4, 2, sharex=True, sharey='row', figsize=(6.5, 5.5))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows):
    for ind, Olim in enumerate(Olims):
        name = common.get_marcher_plot_name(Olim)
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

axes[-1, 0].set_xlabel('Time (s.)')    
axes[-1, 1].set_xlabel('Time (s.)')    

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=3)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.8625)
fig.savefig('time_vs_error_2d.eps')
fig.show()

# size vs error

fig, axes = plt.subplots(4, 2, sharex=True, sharey='row', figsize=(6.5, 5.5))

axes[0, 0].set_title('Relative $\ell_2$ Error')
axes[0, 1].set_title('Relative $\ell_\infty$ Error')

for row, slow in enumerate(Slows):
    for ind, Olim in enumerate(Olims):
        name = common.get_marcher_plot_name(Olim)
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

axes[-1, 0].set_xlabel('$N$')

xticklabels = ['$2^{%d} + 1$' % p for p in Npows]
axes[-1, 1].set_xlabel('$N$')    
axes[-1, 1].set_xticks(N[::2])
axes[-1, 1].set_xticklabels(xticklabels[::2])

handles, labels = axes[-1, -1].get_legend_handles_labels()
    
fig.legend(handles, labels, loc='upper center', ncol=3)
fig.tight_layout()
fig.subplots_adjust(0.05, 0.075, 0.995, 0.8625)
fig.savefig('size_vs_error_2d.eps')
fig.show()
