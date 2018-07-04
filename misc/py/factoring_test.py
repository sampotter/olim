import sys
# sys.path.insert(0, '../../build/Debug')
sys.path.insert(0, '../../build/Release')

import matplotlib.pyplot as plt
plt.ion()

from matplotlib.colors import LogNorm

import numpy as np
import pyeikonal as eik
import speedfuncs

n = 101
h = 2./(n-1)
i0, j0 = int(n/2), int(n/2)
rfac = 1

s = speedfuncs.s5
f = speedfuncs.f5

L = np.linspace(-1, 1, n)
x, y = np.meshgrid(L, L)
S = s(x, y)
F = f(x, y)

r = np.sqrt(x**2 + y**2)
I, J = np.where(r <= rfac)

Olim = eik.Olim8Mid0

m_fac = Olim(S, h)
for i, j in zip(I, J):
    m_fac.set_node_parent(i, j, i0, j0)
m_fac.addBoundaryNode(i0, j0)
m_fac.run()

U_fac = np.array([
    [m_fac.getValue(i, j) for j in range(n)] for i in range(n)])

m = Olim(S, h)
m.addBoundaryNode(i0, j0)
m.run()

U = np.array([
    [m.getValue(i, j) for j in range(n)] for i in range(n)])

eps = np.finfo(np.float64).eps

fig = plt.figure()

fig.add_subplot(3, 2, 1).imshow(U)
fig.add_subplot(3, 2, 2).imshow(U_fac)

D = -np.log(np.maximum(np.abs(U - F), eps))
D_fac = -np.log(np.maximum(np.abs(U_fac - F), eps))
vmin, vmax = min(D.min(), D_fac.min()), max(D.max(), D_fac.max())
norm = LogNorm(vmin=vmin, vmax=vmax)

ax = fig.add_subplot(3, 2, 3)
im = ax.pcolor(D, norm=norm, cmap='PuBu_r')
fig.colorbar(im, ax=ax)

ax = im = fig.add_subplot(3, 2, 4)
im = ax.pcolor(D_fac, norm=norm, cmap='PuBu_r')
fig.colorbar(im, ax=ax)

fig.add_subplot(3, 2, 5).plot(L, U[i0, :])
fig.add_subplot(3, 2, 6).plot(L, U_fac[i0, :])
fig.show()
