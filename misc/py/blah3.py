import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs

norm = np.linalg.norm

plt.ion()
plt.style.use('bmh')

Olim = eik.Olim8Rect
n = 2**10 + 1

h = 2/(n - 1)
i0, j0 = n//2, n//2
l = np.linspace(-1, 1, n)
x, y = np.meshgrid(l, l)

rad = 2*int(np.ceil(np.log2(n)))
dI, dJ = np.meshgrid(np.arange(-rad, rad + 1), np.arange(-rad, rad + 1))
dI = dI.flatten()
dJ = dJ.flatten()

print('number of factored points: %d' % len(I))

u, s = speedfuncs.get_fields(speedfuncs.f0, speedfuncs.s0, x, y)

o = Olim(s, h)
o.addBoundaryNode(i0, j0)
o.run()
U = np.array([[o.getValue(i, j) for i in range(n)] for j in range(n)])

ofac = Olim(s, h)
for di, dj in zip(dI, dJ):
    ofac.set_node_parent(i0 + di, j0 + dj, i0, j0)
ofac.addBoundaryNode(i0, j0)
ofac.run()
Ufac = np.array([[ofac.getValue(i, j) for i in range(n)] for j in range(n)])
print('max difference in circle: %g' %
      np.abs(u[i0 + dI, j0 + dJ] - U[i0 + dI, j0 + dJ]).max())

print('rel l2 err: %g' % (norm(u - U, 'fro')/norm(u, 'fro'),))
print('fac rel l2 err: %g' % (norm(u - Ufac, 'fro')/norm(u, 'fro'),))

fig, ax = plt.subplots(1, 2)
ax[0].grid(False)
im = ax[0].imshow(u - U, cmap=plt.cm.magma, interpolation='none')
fig.colorbar(im, ax=ax[0])
ax[1].grid(False)
im = ax[1].imshow(np.abs(U - Ufac), cmap=plt.cm.magma, interpolation='none')
fig.colorbar(im, ax=ax[1])
plt.show()
