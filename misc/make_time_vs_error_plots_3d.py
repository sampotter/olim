#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys

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
    colors = ['b', 'g', 'r', 'c', 'm']
    if 'basic_marcher_3d' in name: return colors[0]
    elif 'olim6' in name: return colors[1]
    elif 'olim18' in name: return colors[2]
    elif 'olim26' in name: return colors[3]
    elif 'olim3d_hu' in name: return colors[4]

def name_to_label(name):
    if 'basic_marcher_3d' in name: return 'FMM'
    elif 'olim6' in name: return 'OLIM6'
    elif 'olim18' in name: return 'OLIM18'
    elif 'olim26' in name: return 'OLIM26'
    elif 'olim3d_hu' in name: return 'HU OLIM'

def make_loglog_plots(ax, quad, rms_or_inf):
    ax.loglog(data['basic_marcher_3d']['t'],
              data['basic_marcher_3d'][rms_or_inf],
              color=name_to_color('basic_marcher_3d'),
              linewidth=1,
              label=name_to_label('basic_marcher_3d'),
              marker='.')
    for name in data.keys():
        if quad in name:
            ax.loglog(data[name]['t'], data[name][rms_or_inf],
                      linewidth=1, label=name_to_label(name),
                      color=name_to_color(name), marker='.')

if __name__ == '__main__':
    lines = sys.stdin.readlines()
    name_line_indices = [i for i in range(len(lines)) if name_line(lines[i])]
    num_lines = name_line_indices[1] - name_line_indices[0]

    data = dict()
    for i in name_line_indices:
        marcher_name, ns, ts, infs, rmss = \
            parse_csv_section(lines[i:(i + num_lines)])
        data[marcher_name] = {'n': ns, 't': ts, 'inf': infs, 'rms': rmss}

    fig, ax = plt.subplots(3, 2, sharex=True, sharey=True, dpi=300,
                           figsize=(6.5, 6.5))

    make_loglog_plots(ax[0, 0], 'rhr', 'inf')
    make_loglog_plots(ax[1, 0], 'mp0', 'inf')
    make_loglog_plots(ax[2, 0], 'mp1', 'inf')

    make_loglog_plots(ax[0, 1], 'rhr', 'rms')
    make_loglog_plots(ax[1, 1], 'mp0', 'rms')
    make_loglog_plots(ax[2, 1], 'mp1', 'rms')

    ax[0, 0].set_ylabel(r'$F^{(rhr)}$')
    ax[1, 0].set_ylabel('MP0')
    ax[2, 0].set_ylabel('MP1')
    
    ax[0, 0].set_title(r'Relative $\ell_\infty$ Error')
    ax[0, 1].set_title(r'RMS Error')

    # ax[0, 1].set_ylabel(r'Error')
    # ax[0, 0].set_xlabel(r'Time')

    fig.subplots_adjust(
        top=0.95,
        left=0.11,
        right=0.97,
        bottom=0.09)

    fig.legend(loc='lower center', fancybox=False, shadow=False, ncol=5)

    fig.savefig('tmp.pdf')
        
