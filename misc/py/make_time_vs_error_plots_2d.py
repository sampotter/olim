#!/usr/bin/env python3

import glob
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

from itertools import product

def is_csv_line(line):
    return ', ' in line

def name_line(line):
    return not is_csv_line(line)

def parse_csv_line(line):
    n_str, t_str, inf_str, rms_str = line.strip().split(', ')
    return int(n_str), float(t_str), float(inf_str), float(rms_str)

def parse_csv_section(lines):
    marcher_name = lines[0].strip()
    ns, ts, infs, rmss = [], [], [], []
    for n, t, inf, rms in map(parse_csv_line, lines[1:]):
        ns.append(n)
        ts.append(t)
        infs.append(inf)
        rmss.append(rms)
    return marcher_name, ns, ts, infs, rmss

def name_to_color(name):
    return {
        'basic_marcher': 'b',
        'olim4_rhr': 'g',
        'olim4_mp0': 'r',
        'olim4_mp1': 'c',
        'olim8_rhr': 'm',
        'olim8_mp0': 'y',
        'olim8_mp1': 'k'
    }[name]

def name_to_label(name):
    return {
        'basic_marcher': 'FMM',
        'olim4_rhr': 'OLIM4 RHR',
        'olim4_mp0': 'OLIM4 MP0',
        'olim4_mp1': 'OLIM4 MP1',
        'olim8_rhr': 'OLIM8 RHR',
        'olim8_mp0': 'OLIM8 MP0',
        'olim8_mp1': 'OLIM8 MP1'
    }[name]

def make_loglog_plots(ax, rms_or_inf):
    for name in data.keys():
        ax.loglog(data[name]['t'], data[name][rms_or_inf],
                  linewidth=1, label=name_to_label(name),
                  color=name_to_color(name), marker='.')

if __name__ == '__main__':

    data = dict()

    paths = glob.glob('../data/plot_data_2d.*.txt')

    def get_sname(path):
        return path.split('/')[-1].split('.')[1]

    for path in paths:
        sname = get_sname(path)

        with open(path) as f:
            lines = f.readlines()

            name_line_indices = \
                [i for i in range(len(lines)) if name_line(lines[i])]
            num_lines = name_line_indices[1] - name_line_indices[0]

        data[sname] = dict()
        for i in name_line_indices:
            marcher_name, ns, ts, infs, rmss = \
                parse_csv_section(lines[i:(i + num_lines)])
            data[sname][marcher_name] = \
                {'n': ns, 't': ts, 'inf': infs, 'rms': rmss}

    snames = [get_sname(p) for p in paths]
    if 's2' in snames: snames.remove('s2')

    olims = list(data[snames[0]].keys())

    nrows, ncols = 4, 4

    fig, axes = plt.subplots(nrows, ncols, sharex=True, sharey=True, dpi=300,
                           figsize=(6.5, 8))

    for i, j in product(range(int(nrows/2)), range(ncols)):

        sname = snames[ncols*i + j]

        ax1 = axes[i, j]
        ax2 = axes[int(nrows/2) + i, j]

        for olim in olims:
            d = data[sname][olim]
            ax1.loglog(d['t'], d['inf'], label=olim.replace('_', ' '),
                       linewidth=1, marker='|', markersize=3.5)
            ax2.loglog(d['t'], d['rms'], label=olim.replace('_', ' '),
                       linewidth=1, marker='|', markersize=3.5)

        ax1.text(0.05, 0.05, r'\texttt{%s}' % sname, transform=ax1.transAxes)
        ax2.text(0.05, 0.05, r'\texttt{%s}' % sname, transform=ax2.transAxes)

    fig.tight_layout()
    # fig.show()
    fig.savefig('../data/plots_2d.eps')
        
