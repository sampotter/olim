#!/usr/bin/env python3

import fileinput
import numpy as np

def rhs2int(s):
    return int(s.split('=')[1])

def rhs2ratio(s):
    return tuple(map(int, s.split('=')[1].split('/')))

if __name__ == '__main__':
    file_in = fileinput.input()

    first_line = file_in.readline()
    depth, width, height = map(rhs2int, first_line.split(','))

    visits = np.empty((depth, width, height), dtype=np.int32)
    lines = np.empty_like(visits)
    tris = np.empty_like(visits)
    tetras = np.empty_like(visits)

    for line in file_in:
        lhs, rhs = line.split(':')
        i, j, k = map(int, lhs.split(','))
        tmp = rhs.split(',')
        visits[i, j, k] = rhs2int(tmp[0])
        lines[i, j, k] = rhs2int(tmp[1])
        tris[i, j, k] = rhs2ratio(tmp[2])[1]
        tetras[i, j, k] = rhs2ratio(tmp[3])[1]

    print('<line> = %g, <tri> = %g, <tetra> = %g (<visits> = %g)' % (
        lines.mean(),
        tris.mean(),
        tetras.mean(),
        visits.mean()))
