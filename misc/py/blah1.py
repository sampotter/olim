import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs3d as s3d

Olim = eik.Olim26Mid1
r_fac = 0.1
N = 2**np.arange(3, 9) + 1

for n in N:

    h = 2/(n-1)
    i0 = n//2

    l = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(l, l, l)
    r = np.sqrt(x**2 + y**2 + z**2)

    s = s3d.s1(x, y, z)
    u = s3d.f1(x, y, z)

    o = Olim(s, h)
    o.addBoundaryNode(i0, i0, i0)
    for i, j, k in zip(*np.where(r <= r_fac)):
        o.set_node_parent(i, j, k, i0, i0, i0)
    o.run()

    U = np.array([[[o.getValue(i, j, k) for j in range(n)]
                   for i in range(n)]
                  for k in range(n)])

    err2 = np.linalg.norm((u - U).flatten())/np.linalg.norm(u.flatten())

    print(err2)
