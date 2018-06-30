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
    lines = sys.stdin.readlines()
    name_line_indices = [i for i in range(len(lines)) if name_line(lines[i])]
    num_lines = name_line_indices[1] - name_line_indices[0]

    data = dict()
    for i in name_line_indices:
        marcher_name, ns, ts, infs, rmss = \
            parse_csv_section(lines[i:(i + num_lines)])
        data[marcher_name] = {'n': ns, 't': ts, 'inf': infs, 'rms': rmss}

    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, dpi=300,
                           figsize=(6.5, 6.5))

    make_loglog_plots(ax[0], 'inf')
    make_loglog_plots(ax[1], 'rms')

    ax[0].set_title(r'Relative $\ell_\infty$ Error')
    ax[1].set_title(r'RMS Error')

    # ax[0, 1].set_ylabel(r'Error')
    # ax[0, 0].set_xlabel(r'Time')

    fig.subplots_adjust(
        top=0.95,
        left=0.08,
        right=0.98,
        bottom=0.16)

    fig.legend(loc='lower center', fancybox=False, shadow=False, ncol=3)

    fig.savefig('tmp.pdf')
        
