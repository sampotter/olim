import sys
if '../../build/Release' not in sys.path:
    sys.path.insert(0, '../../build/Release')

import matplotlib.pyplot as plt
import numpy as np
import pyeikonal as eik
import speedfuncs3d

plt.ion()
plt.style.use('classic')

# Olim = eik.BasicMarcher3D
Olim = eik.Olim18Rect
N = 2**np.arange(3, 8) + 1

for n in N:

    h = 2/(n-1)
    i0 = n//2

    l = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(l, l, l)

    u, s = speedfuncs3d.get_fields(speedfuncs3d.f4, speedfuncs3d.s4, x, y, z)

    o = Olim(s, h)
    o.addBoundaryNode(i0, i0, i0)
    o.run()

    U = np.array([[[o.getValue(i, j, k) for i in range(n)] for j in range(n)]
                  for k in range(n)])

    err2 = np.linalg.norm((u - U).flatten())/np.linalg.norm(u.flatten())

    print(err2)

# fig, axes = plt.subplots(3, 2)
# axes[0, 0].contourf(x[:, :, i0], y[:, :, i0], u[:, :, i0])
# axes[0, 1].contourf(x[:, :, i0], y[:, :, i0], U[:, :, i0])
# axes[1, 0].contourf(x[:, :, i0], y[:, :, i0], u[:, i0, :])
# axes[1, 1].contourf(x[:, :, i0], y[:, :, i0], U[:, i0, :])
# axes[2, 0].contourf(x[:, :, i0], y[:, :, i0], u[i0, :, :])
# axes[2, 1].contourf(x[:, :, i0], y[:, :, i0], U[i0, :, :])
# plt.colorbar()
# plt.show()

sl = -1
fig, axes = plt.subplots(1, 3, sharey=True)
cs = axes[0].contourf(x[:,:,i0], y[:,:,i0], (u - U)[:,:,sl], 10)
plt.colorbar(cs, ax=axes[0])
cs = axes[1].contourf(x[:,:,i0], y[:,:,i0], u[:,:,sl], 10)
plt.colorbar(cs, ax=axes[1])
cs = axes[2].contourf(x[:,:,i0], y[:,:,i0], s[:,:,sl], 10)
plt.colorbar(cs, ax=axes[2])
plt.tight_layout()
plt.show()
