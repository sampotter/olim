#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--min_3d_power', type=int, default=3)
parser.add_argument('--max_3d_power', type=int, default=15)
parser.add_argument('--build_type', type=str, default='Release')
args = parser.parse_args()

################################################################################
# Preliminaries

import glob
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'serif': ['Computer Modern']})

plt.style.use('bmh')

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

################################################################################
# Main

olim_strs = ['olim6', 'olim18', 'olim26', 'olim3d_hu']
P = np.arange(args.min_3d_power, args.max_3d_power + 1)
N = 2**P + 1

visits = {k: [] for k in olim_strs}
line = {k: [] for k in olim_strs}
tri = {k: [] for k in olim_strs}
tetra = {k: [] for k in olim_strs}

# Collect statistics

for olim_str, n in it.product(olim_strs, N):
    print('%s, n = %d' % (olim_str, n))
    subprocess.run([
        '../build/Stats/gen_stats', olim_str + '_rhr', str(n), 'stats.bin'])

    S = np.fromfile('stats.bin', np.int32).reshape(n**3, 7)

    visits[olim_str].append(S[:, 3].mean())
    line[olim_str].append(S[:, 4].mean())
    tri[olim_str].append(S[:, 5].mean())
    tetra[olim_str].append(S[:, 6].mean())

################################################################################
# Plotting
    
fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(6.5, 4))
ax = axes.flatten()

for i, olim_str in enumerate(olim_strs):
    ax[i].semilogx(N, line[olim_str], linewidth='1', color=colors[i],
                   linestyle='dotted', marker='|', markersize=3.5,
                   label=r'$d = 0$')
    ax[i].semilogx(N, tri[olim_str], linewidth='1', color=colors[i],
                   linestyle='dashed', marker='|', markersize=3.5,
                   label=r'$d = 1$')
    ax[i].semilogx(N, tetra[olim_str], linewidth='1', color=colors[i],
                   linestyle='solid', marker='|', markersize=3.5,
                   label=r'$d = 2$')
    loc = 'upper left' if i == 0 else 'lower left'
    ax[i].legend(loc=loc, ncol=1, prop={'size': 8})
    ax[i].minorticks_off()
    ax[i].text(0.9, 0.9 if i == 0 else 0.1,
               r'$\texttt{%s}$' % olim_str.replace('_', '\_'),
               transform=ax[i].transAxes, horizontalalignment='center',
               verticalalignment='center')

for j in range(2):
    axes[1, j].set_xticks(2**P[::3] + 1)
    axes[1, j].set_xticklabels(['$2^{%d} + 1$' % p for p in P[::3]])
    axes[1, j].set_xlabel(r'$N$')
    
fig.tight_layout()
fig.show()

fig.savefig('stats.eps')
