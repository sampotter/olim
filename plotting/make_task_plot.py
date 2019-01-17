#!/usr/bin/env python3

################################################################################
# parse arguments first

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument('--min_2d_power', type=int, default=3)
    # parser.add_argument('--max_2d_power', type=int, default=15)
    # parser.add_argument('--build_type', type=str, default='Release')
    parser.add_argument('path', type=str)
    args = parser.parse_args()

################################################################################
# preliminaries

import fileinput

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


def insert_into_data_frame(d, df, N):
    update_keys = {'update_impl'}
    control_keys = {'visit_neighbors_impl', 'run'}
    heap_keys = {'adjust_heap_entry', 'insert_into_heap'}
    other_keys = {'unknown', 'insert', 'marcher_3d'}
    for k in d:
        if k in update_keys:
            df['update'][N] += d[k]
        elif k in control_keys:
            df['control'][N] += d[k]
        elif k in heap_keys:
            df['heap'][N] += d[k]
        elif k in other_keys:
            df['other'][N] += d[k]
        else:
            raise Exception('missing key: %s' % k)
