#!/usr/bin/env python3

import fileinput
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from matplotlib.path import Path

def line_seg(p, q):
    return Path([p, q], [Path.MOVETO, Path.LINETO])

def convcomb(p0, p1, lam):
    return (1 - lam)*p0 + lam*p1

class node(object):
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.T = None
        self.lam = None
        self.i0 = None
        self.j0 = None
        self.i1 = None
        self.j1 = None

    def get_type(self):
        if self.i1:
            return 'TWO_POINT'
        elif self.i0:
            return 'ONE_POINT'
        else:
            return 'BOUNDARY'

    def get_path(self):
        type = self.get_type()
        if type == 'TWO_POINT':
            ilam = convcomb(self.i0, self.i1, self.lam)
            jlam = convcomb(self.j0, self.j1, self.lam)
            return Path.make_compound_path(
                line_seg((self.i0, self.j0), (self.i1, self.j1)),
                line_seg((self.i, self.j), (ilam, jlam)),
                Path.circle((self.i, self.j), 0.1),
                Path.circle((ilam, jlam), 0.05))
        elif type == 'ONE_POINT':
            return Path.make_compound_path(
                line_seg((self.i, self.j), (self.i0, self.j0)),
                Path.circle((self.i, self.j), 0.1))
        else:
            return Path.circle((self.i, self.j), 0.1)

    def __str__(self):
        type = self.get_type()
        if type == 'TWO_POINT':
            return 'node(%d, %d, %g, %g, %d, %d, %d, %d)' % (
                self.i, self.j, self.T, self.lam, self.i0, self.j0, self.i1,
                self.j1)
        elif type == 'ONE_POINT':
            return 'node(%d, %d, %d, %d, %d)' % (
                self.i, self.j, self.T, self.i0, self.j0)
        else:
            return 'node(%d, %d, %d)' % (self.i, self.j, self.T)

def parse_node_from_line(str):
    start, end = str.split(':')
    i, j = map(int, start.split(','))
    n = node(i, j)
    for k, v in map(lambda s: s.split('='), end.split(',')):
        k = k.strip()
        if k == 'T': n.T = float(v)
        elif k == 'lam': n.lam = float(v)
        elif k == 'i0': n.i0 = int(v)
        elif k == 'j0': n.j0 = int(v)
        elif k == 'i1': n.i1 = int(v)
        elif k == 'j1': n.j1 = int(v)
        else: raise Exception('key error: got %s' % k)
    return n

def get_bounds(nodes):
    imin, imax = float('+inf'), float('-inf')
    jmin, jmax = float('+inf'), float('-inf')
    for n in nodes:
        imin = min(imin, n.i)
        imax = max(imax, n.i)
        jmin = min(jmin, n.j)
        jmax = max(jmax, n.j)
    return imin, imax, jmin, jmax

def get_grid_path(imin, imax, jmin, max):
    hlines, vlines = [], []
    for j in range(jmin, jmax + 1):
        hlines.append(line_seg((imin, j), (imax, j)))
    for i in range(imin, imax + 1):
        vlines.append(line_seg((i, jmin), (i, jmax)))
    return Path.make_compound_path(*(hlines + vlines))

if __name__ == '__main__':
    nodes = []
    for line in fileinput.input():
        if line.find(':') != -1:
            nodes.append(parse_node_from_line(line))
    imin, imax, jmin, jmax = get_bounds(nodes)

    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    ax.set_xlim(jmin - 1, jmax + 1)
    ax.set_ylim(imin - 1, imax + 1)

    path = get_grid_path(imin, imax, jmin, jmax)
    patch = patches.PathPatch(path, facecolor='gray')
    ax.add_patch(patch)

    for n in nodes:
        path = n.get_path()
        patch = patches.PathPatch(path, facecolor='black', lw=2.5)
        ax.add_patch(patch)
    
    plt.show()

# Local Variables:
# indent-tabs-mode: nil
# End:
