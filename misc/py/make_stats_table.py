#!/usr/bin/env python3

import argparse
import glob
import numpy as np
import os.path

np.seterr(invalid='ignore') # don't care if we divide by zero in this script

DATA_DIR = '../data'
STATS_FILES = glob.glob(os.path.join(DATA_DIR, 'olim*.txt'))

def get_problem_shape(header):
    return tuple(int(s.split('=')[1]) for s in header.split(','))

def process_line(line):
    ijk_str, line = line.split(':')
    i, j, k = map(int, ijk_str.split(','))
    visits_str, line_str, tri_str, tetra_str = line.split(',')
    return np.array([
        int(visits_str.split('=')[1]),
        int(line_str.split('=')[1]), 
        *tuple(map(int, tri_str.split('=')[1].split('/'))),
        *tuple(map(int, tetra_str.split('=')[1].split('/')))])

def init_stats_arrays(n):
    return np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n)), \
        np.empty((n, n, n)), np.empty((n, n, n)), np.empty((n, n, n))

def load_stats_file(path):
    with open(path) as f:
        lines = f.readlines()
    header, lines = lines[0], lines[1:]
    d, w, h = get_problem_shape(header)
    return np.array(list(map(process_line, lines))).reshape(d, w, h, 6)

def stats_file_path_to_size(s):
    return int(s.split('/')[-1].split('.')[1])

def stats_file_path_to_olim_and_size(s):
    olim, size_str = s.split('/')[-1].split('.')[:2]
    size = int(size_str)
    return olim, size

def build_table(max_n=np.inf, quiet=False):
    table = dict()
    for path in STATS_FILES:
        olim, n = stats_file_path_to_olim_and_size(path)
        if n > max_n or 'rhr' not in olim:
            continue
        print('- %s, n = %d' % (olim, n))
        if olim not in table:
            table[olim] = dict()
        if n not in table[olim]:
            table[olim][n] = dict()
        def get_ratio_mean(i):
            return np.nanmean(stats[:, :, :, i]/stats[:, :, :, 0])
        stats = load_stats_file(path)
        table[olim][n]['avg_visits'] = stats[:,:,:,0].mean()
        table[olim][n]['avg_line'] = get_ratio_mean(1)
        table[olim][n]['avg_tri_nd'] = get_ratio_mean(2)
        table[olim][n]['avg_tri_tot'] = get_ratio_mean(3)
        table[olim][n]['avg_tetra_nd'] = get_ratio_mean(4)
        table[olim][n]['avg_tetra_tot'] = get_ratio_mean(5)
    return table

def make_tabular_lines(table):
    olims = list(table.keys())
    ns = sorted(table[olims[0]].keys())
    keys = list(table[olims[0]][ns[0]].keys())

    lines = []
    lines.append(r'\begin{tabular}{c|r|r|r|rr|rr}')
    lines.append('&' + ' & '.join([
        r'\multirow{2}{*}{\centering $n$}',
        r'\multirow{2}{*}{Avg. Visits}',
        r'\multirow{2}{*}{$d = 0$}',
        r'\multicolumn{2}{c}{$d = 1$}',
        r'\multicolumn{2}{c}{$d = 2$}']) + r' \\')
    lines.append(
        ' & '.join(
            ['& & & '] + 
            ([r'\multicolumn{1}{c}{Nondeg.}', r'\multicolumn{1}{c}{Total}'] * 2)) +
        r' \\')
    for olim in olims:
        lines.append(r'\midrule')
        for i, n in enumerate(ns):
            line = (r'& %d & ' % n) + ' & '.join([
                '%0.4f' % x for x in [table[olim][n][key] for key in keys]]) + \
                r' \\'
            if i == 0:
                line = (r'\multirow{%d}{*}{\texttt{%s}} ' % (
                    len(ns), olim.replace('_', r'\_'))) + line
            lines.append(line)
    lines.append(r'\midrule')
    lines.append(r'\end{tabular}')

    return lines

def output_table(table, output_path):
    with open(output_path, 'w') as f:
        for line in make_tabular_lines(table):
            print(line, file=f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output_path', type=str)

    args = parser.parse_args()
    output_path = args.output_path

    table = build_table()
    output_table(table, output_path)
