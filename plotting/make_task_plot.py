#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse
import os

build_dir = '../build/RelWithDebInfo'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--min_3d_power', type=int, default=3)
    parser.add_argument('--max_3d_power', type=int, default=15)
    parser.add_argument('--build_type', type=str, default='RelWithDebInfo')
    parser.add_argument('--tmp_dir', type=str, default='./tmp')
    args = parser.parse_args()

    build_dir = os.path.join('..', 'build', args.build_type)

################################################################################
# preliminaries

import sys
sys.path.insert(0, '../build/%s' % args.build_type)
sys.path.insert(0, '../misc/py')

import fileinput
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyolim as olim
import subprocess

from matplotlib.lines import Line2D

from cycler import cycler
from matplotlib import rc

rc('text', usetex=True)
rc('font', **{
    'family': 'serif',
    'serif': ['Computer Modern'],
    'size': 8
})

plt.ion()
plt.style.use('bmh')

################################################################################
# Collect data

def parse_comma_sep_int(s):
    return int(s.replace(',', ''))

def parse_func_name(s):
    l, r = s.find('::'), s.rfind('(')
    if l < 0 or s.find('scratch') < 0:
        return 'unknown'
    else:
        l += 2
        s = s[l:r]
        r = s.find('<')
        if r >= 0:
            s = s[:r]
        r = s.find('(')
        if r >= 0:
            s = s[:r]
        if s.find('__') >= 0:
            s = 'unknown'
        return s

def add_or_increment(d, k, v):
    if k not in d:
        d[k] = 0
    d[k] += v

def parse_ir_dict(path):

    f = fileinput.input(path)

    func2ir = dict()

    while True:
        s = f.readline()
        if s is None or s.strip() == 'Ir':
            break

    f.readline()

    s = f.readline().split()[0]
    total_ir = parse_comma_sep_int(s)

    for _ in range(4): f.readline()

    s = f.readline().strip()
    while len(s) > 0:
        func = parse_func_name(s)
        ir = parse_comma_sep_int(s.split()[0])
        add_or_increment(func2ir, func, ir)
        s = f.readline().strip()

    ir_diff = total_ir - sum(func2ir.values()) 
    assert(ir_diff >= 0)
    if ir_diff > 0:
        add_or_increment(func2ir, 'unknown', ir_diff)
    assert(sum(func2ir.values()) == total_ir)

    f.close()

    return func2ir

def make_cg_annotation(Olim, name, n, path):
    subprocess.run([
        'valgrind',
        '--tool=callgrind',
        '--callgrind-out-file=callgrind.out',
        os.path.join(build_dir, 'scratch'),
        name,
        str(n)])
    with open(path, 'w') as f:
        subprocess.run(['callgrind_annotate', 'callgrind.out'], stdout=f)

Olims = [olim.BasicMarcher3D, olim.Olim6Rect, olim.Olim26Mid0, olim.Olim26Mid1,
         olim.Olim3dHuMid0, olim.Olim3dHuMid1]

names = ['basic_marcher_3d', 'olim6_rhr', 'olim26_mp0', 'olim26_mp1',
         'olim3d_hu_mp0', 'olim3d_hu_mp1']

P = np.arange(args.min_3d_power, args.max_3d_power + 1)
N = 2**P + 1

assert(len(Olims) == len(names))

os.makedirs(args.tmp_dir, exist_ok=True)
for (Olim, name), n in it.product(zip(Olims, names), N):
    print('\n%s (n = %d)\n' % (name, n))
    path = os.path.join(args.tmp_dir, '%s_%d.txt' % (name, n))
    make_cg_annotation(Olim, name, n, path)

def initialize_or_increment(df, col, row, val):
    if np.isnan(df[col][row]):
        df[col][row] = val
    else:
        df[col][row] += val

def insert_into_data_frame(d, df, N):
    update_keys = {'update_crtp', 'update_impl', 'tetra', 'tri', 'tri_bv',
                   'operator', 'should_skip'}
    heap_keys = {'adjust_heap_entry', 'insert_into_heap', 'swim',
                 'sink', 'get_next_node'}
    logic_keys = {'visit_neighbors_impl', 'run', 'unknown', 'insert',
                  'marcher_3d', 'init', 'init_crtp', '_Function_handler',
                  'conditional', 'pair'}
    for k in d:
        if k in update_keys:
            initialize_or_increment(df, 'update', N, d[k])
        elif k in logic_keys:
            initialize_or_increment(df, 'logic', N, d[k])
        elif k in heap_keys:
            initialize_or_increment(df, 'heap', N, d[k])
        else:
            raise Exception('missing key: %s' % k)

task_dfs = {Olim: pd.DataFrame(index=N,
                               columns=['update', 'heap', 'logic'])
            for Olim in Olims}
for (Olim, name), n in it.product(zip(Olims, names), N):
    print('%s (n = %d)' % (name, n))
    path = os.path.join(args.tmp_dir, '%s_%d.txt' % (name, n))
    d = parse_ir_dict(path)
    insert_into_data_frame(d, task_dfs[Olim], n)

for Olim in Olims:
    df = task_dfs[Olim]
    task_dfs[Olim] = df.divide(df.sum(axis=1), axis=0)
    
################################################################################
# Plot task percentages

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
cmap = [0, 1, 4, 3]
colors = [colors[i] for i in cmap]

style = {
    'linewidth': 1,
    'marker': '|',
    'markersize': 3.5
}

def plot_columns(ax, df, linestyle):
    ax.semilogx(N, df['update'], color=colors[0], linestyle=linestyle, **style)
    ax.semilogx(N, df['heap'], color=colors[1], linestyle=linestyle, **style)
    ax.semilogx(N, df['logic'], color=colors[2], linestyle=linestyle, **style)

fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.5), sharex=True, sharey=True)

ax = axes[0]
plot_columns(ax, task_dfs[Olims[0]], '-')
plot_columns(ax, task_dfs[Olims[1]], '--')
ax.set_ylim(-0.05, 1)
ax.minorticks_off()
ax.set_xticks(N[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in P[::2]])
ax.set_xlabel('$N$')
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels([r'0\%', r'25\%', r'50\%', r'75\%', r'100\%'])

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::3], [r'FMM', r'\texttt{olim6\_rhr}'])

ax = axes[1]
plot_columns(ax, task_dfs[Olims[2]], '-')
plot_columns(ax, task_dfs[Olims[3]], '--')
ax.set_ylim(-0.05, 1)
ax.minorticks_off()
ax.set_xticks(N[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in P[::2]])
ax.set_xlabel('$N$')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::3], [r'\texttt{olim26\_mp0}', r'\texttt{olim26\_mp1}'])

ax = axes[2]
plot_columns(ax, task_dfs[Olims[4]], '-')
plot_columns(ax, task_dfs[Olims[5]], '--')
ax.set_ylim(-0.05, 1)
ax.minorticks_off()
ax.set_xticks(N[::2])
ax.set_xticklabels(['$2^{%d} + 1$' % p for p in P[::2]])
ax.set_xlabel('$N$')

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::3], [r'\texttt{olim3d\_mp0}', r'\texttt{olim3d\_mp1}'])

fig.legend(handles[:3], labels[:3], ncol=3, loc='upper center')

fig.tight_layout()
fig.subplots_adjust(0.05, 0.13, 0.995, 0.92, 0.12, 0.20)
fig.show()

fig.savefig('tasks.eps')

# colors = ['black', 'red', 'green']
# lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='--') for c in colors]
# labels = ['black data', 'red data', 'green data']
# plt.legend(lines, labels)
# plt.show()
