#!/usr/bin/env python3

import glob
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import subprocess

from matplotlib import rc

rc('text', usetex=True)
rc('font', **{'serif': ['Computer Modern']})

plt.style.use('bmh')

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

olim_strs = ['olim6', 'olim18', 'olim26', 'olim3d_hu']
N = 2**np.arange(3, 8) + 1

if __name__ == '__main__':

    visits = {k: [] for k in olim_strs}
    line = {k: [] for k in olim_strs}
    tri = {k: [] for k in olim_strs}
    tetra = {k: [] for k in olim_strs}

    for olim_str, n in it.product(olim_strs, N):
        print('%s, n = %d' % (olim_str, n))
        subprocess.run([
            '../build/Stats/gen_stats', olim_str + '_rhr', str(n), 'stats.bin'])
    
        S = np.fromfile('stats.bin', np.int32).reshape(n**3, 7)

        visits[olim_str].append(S[:, 3].mean())
        line[olim_str].append(S[:, 4].mean())
        tri[olim_str].append(S[:, 5].mean())
        tetra[olim_str].append(S[:, 6].mean())

    plt.figure(figsize=(6, 4))
    for i, olim_str in enumerate(olim_strs):
        # plt.semilogx(N, visits[olim_str], label=olim_str)
        plt.semilogx(N, line[olim_str], linewidth='2',
                     color=colors[i], linestyle='dotted')
        plt.semilogx(N, tri[olim_str], linewidth='2',
                     color=colors[i], linestyle='dashed')
        plt.semilogx(N, tetra[olim_str], linewidth='2',
                     color=colors[i], linestyle='solid',
                     label=olim_str.replace('_', '\_'))
    plt.legend()
    plt.xlabel('$N$')
    plt.show()
