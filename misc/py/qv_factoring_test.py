import sys; sys.path.insert(0, '../../build/Release')

import common
import common3d
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyeikonal as eik

from matplotlib.colors import LogNorm
from numpy.linalg import norm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.style.use('bmh')

################################################################################
# parameters

R_fac = 0.25
N = 2**np.arange(3, 8) + 1
N3D = 2**np.arange(3, 6) + 1
vx, vy, vz = 5, 13, 20

x_fac_1, y_fac_1, z_fac_1 = 0.0, 0.0, 0.0
x_fac_2, y_fac_2, z_fac_2 = 0.8, 0.0, 0.0

################################################################################
# common

def rms(u, U):
    return np.sqrt((np.mean((u - U).flatten())**2))

def linf(u, U):
    return norm((u - U).flatten(), np.inf)/max(
        norm(u.flatten(), np.inf), norm(U.flatten(), np.inf))

################################################################################
# 2D

s = lambda x, y: 1/(2 + vx*x + vy*y)

def make_u(x_fac, y_fac, vx, vy, s):
    return lambda x, y: \
        (1/np.sqrt(vx**2 + vy**2)) * \
        np.arccosh(
            1 +
            s(x_fac, y_fac)*s(x, y)*(vx**2 + vy**2)*(
                (x - x_fac)**2 + (y - y_fac)**2)/2)

u_1 = make_u(x_fac_1, y_fac_1, vx, vy, s)
u_2 = make_u(x_fac_2, y_fac_2, vx, vy, s)

u = lambda x, y: np.minimum(u_1(x, y), u_2(x, y))

E = dict()
E_fac = dict()

for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    
    print(common.get_marcher_name(Olim))

    E[Olim] = {'rms': np.zeros(len(N)), 'inf': np.zeros(len(N))}
    E_fac[Olim] = {'rms': np.zeros(len(N)), 'inf': np.zeros(len(N))}

    for k, n in enumerate(N):
        print('- n = %d (%d/%d)' % (n, k + 1, len(N)))

        L = np.linspace(0, 1, n)
        X, Y = np.meshgrid(L, L)
        u_ = u(X, Y)
        S = s(X, Y)

        h = 1/(n - 1)
        i_1, j_1 = int(y_fac_1/h), int(x_fac_1/h)
        i_2, j_2 = int(y_fac_2/h), int(x_fac_2/h)

        m = Olim(S, h)
        m.addBoundaryNode(i_1, j_1)
        m.addBoundaryNode(i_2, j_2)
        m.run()
        U = np.array([[m.getValue(i, j) for j in range(n)] for i in range(n)])

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2)
        for i, j in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_parent(i, j, i_1, j_1)
        m_fac.addBoundaryNode(i_1, j_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2)
        for i, j in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_parent(i, j, i_2, j_2)
        m_fac.addBoundaryNode(i_2, j_2)

        m_fac.run()
        U_fac = np.array(
            [[m_fac.getValue(i, j) for j in range(n)] for i in range(n)])

        E[Olim]['rms'][k] = rms(u_, U)
        E[Olim]['inf'][k] = linf(u_, U)

        E_fac[Olim]['rms'][k] = rms(u_, U_fac)
        E_fac[Olim]['inf'][k] = linf(u_, U_fac)

# plotting

fig, axes = plt.subplots(
    2, 2, sharex=True, figsize=(6.5, 4.5), dpi=100)

axes[0, 0].set_ylabel(r'Unfactored')
axes[1, 0].set_ylabel(r'Factored')
axes[0, 0].set_title(r'RMS Error')
axes[0, 1].set_title(r'Relative $\ell_\infty$ Error')

ax = axes[0, 0]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E[Olim]['rms'],
        label=common.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

ax = axes[0, 1]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E[Olim]['inf'],
        label=common.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

ax = axes[1, 0]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E_fac[Olim]['rms'],
        label=common.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

ax = axes[1, 1]
for Olim in common.marchers:
    if Olim == eik.BasicMarcher:
        continue
    ax.loglog(
        N, E_fac[Olim]['inf'],
        label=common.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
ax.minorticks_off()
ax.set_xticks(N)
ax.set_xticklabels(map(str, N))

fig.legend()

fig.tight_layout()
fig.show()
# fig.savefig('../data/qv_plots_2d.eps')

################################################################################
# 3D

s3d = lambda x, y, z: 1/(2 + vx*x + vy*y + vz*z)

def make_u3d(x_fac, y_fac, z_fac, vx, vy, vz, s):
    return lambda x, y, z: \
        (1/np.sqrt(vx**2 + vy**2 + vz**2)) * \
        np.arccosh(
            1 +
            s3d(x_fac, y_fac, z_fac)*s3d(x, y, z)*(vx**2 + vy**2 + vz**2)*
            ((x - x_fac)**2 + (y - y_fac)**2 + (z - z_fac)**2)/2)

u3d_1 = make_u3d(x_fac_1, y_fac_1, z_fac_1, vx, vy, vz, s)
u3d_2 = make_u3d(x_fac_2, y_fac_2, z_fac_2, vx, vy, vz, s)

u3d = lambda x, y, z: np.minimum(u3d_1(x, y, z), u3d_2(x, y, z))

E = dict()
E_fac = dict()

# Olims = common3d.marchers.copy()
# Olims.remove(eik.BasicMarcher3D)
# Olims.remove(eik.Olim3dHuMid0)
# Olims.remove(eik.Olim3dHuMid1)
# Olims.remove(eik.Olim3dHuRect)

Olims = [eik.BasicMarcher3D]

for Olim in Olims:
    print(common3d.get_marcher_name(Olim))

    E[Olim] = {'rms': np.zeros(len(N3D)), 'inf': np.zeros(len(N3D))}
    E_fac[Olim] = {'rms': np.zeros(len(N3D)), 'inf': np.zeros(len(N3D))}

    for a, n in enumerate(N3D):
        print('- n = %d (%d/%d)' % (n, a + 1, len(N3D)))

        L = np.linspace(0, 1, n)
        X, Y, Z = np.meshgrid(L, L, L)
        u_ = u3d(X, Y, Z)
        S = s3d(X, Y, Z)

        h = 1/(n - 1)
        i_1, j_1, k_1 = int(y_fac_1/h), int(x_fac_1/h), int(z_fac_1/h)
        i_2, j_2, k_2 = int(y_fac_2/h), int(x_fac_2/h), int(z_fac_2/h)
        # i_3, j_3, k_3 = int(y_fac_3/h), int(x_fac_3/h), int(z_fac_3/h)

        m = Olim(S, h)
        m.addBoundaryNode(i_1, j_1, k_1)
        m.addBoundaryNode(i_2, j_2, k_2)
        # m.addBoundaryNode(i_3, j_3, k_3)
        m.run()
        U = np.array([[[m.getValue(i, j, k) for i in range(n)]
                       for j in range(n)]
                      for k in range(n)])

        m_fac = Olim(S, h)

        R_1 = np.sqrt((x_fac_1 - X)**2 + (y_fac_1 - Y)**2 + (z_fac_1 - Z)**2)
        for i, j, k in zip(*np.where(R_1 <= R_fac)):
            m_fac.set_node_parent(i, j, k, i_1, j_1, k_1)
        m_fac.addBoundaryNode(i_1, j_1, k_1)

        R_2 = np.sqrt((x_fac_2 - X)**2 + (y_fac_2 - Y)**2 + (z_fac_2 - Z)**2)
        for i, j, k in zip(*np.where(R_2 <= R_fac)):
            m_fac.set_node_parent(i, j, k, i_2, j_2, k_2)
        m_fac.addBoundaryNode(i_2, j_2, k_2)

        # R_3 = np.sqrt((x_fac_3 - X)**2 + (y_fac_3 - Y)**2 + (z_fac_3 - Z)**2)
        # for i, j, k in zip(*np.where(R_3 <= R_fac)):
        #     m_fac.set_node_parent(i, j, k, i_3, j_3, k_3)
        # m_fac.addBoundaryNode(i_3, j_3, k_3)

        m_fac.run()
        U_fac = np.array([[[m_fac.getValue(i, j, k) for i in range(n)]
                           for j in range(n)]
                          for k in range(n)])

        E[Olim]['rms'][a] = rms(u_, U)
        E[Olim]['inf'][a] = linf(u_, U)

        E_fac[Olim]['rms'][a] = rms(u_, U_fac)
        E_fac[Olim]['inf'][a] = linf(u_, U_fac)

# debugging

X_ = X[:, :, 0]
Y_ = Y[:, :, 0]
fig, ax = plt.subplots(3, 2, sharex=True, sharey=True)
ax[0, 0].set_title('z front')
fig.colorbar(ax[0, 0].contourf(X_, Y_, (U - u_)[:, :, 0]), ax=ax[0, 0])
ax[0, 1].set_title('z back')
fig.colorbar(ax[0, 1].contourf(X_, Y_, (U - u_)[:, :, -1]), ax=ax[0, 1])
ax[1, 0].set_title('y front')
fig.colorbar(ax[1, 0].contourf(X_, Y_, (U - u_)[:, 0, :]), ax=ax[1, 0])
ax[1, 1].set_title('y back')
fig.colorbar(ax[1, 1].contourf(X_, Y_, (U - u_)[:, -1, :]), ax=ax[1, 1])
ax[2, 0].set_title('x front')
fig.colorbar(ax[2, 0].contourf(X_, Y_, (U - u_)[0, :, :]), ax=ax[2, 0])
ax[2, 1].set_title('x back')
fig.colorbar(ax[2, 1].contourf(X_, Y_, (U - u_)[-1, :, :]), ax=ax[2, 1])
fig.show()

# X_ = X[:, :, 0]
# Y_ = Y[:, :, 0]
# fig, ax = plt.subplots(3, 2, sharex=True, sharey=True)
# ax[0, 0].contourf(X_, Y_, u_[:, :, 0])
# ax[0, 1].contourf(X_, Y_, U[:, :, 0])
# ax[1, 0].contourf(X_, Y_, u_[:, 0, :])
# ax[1, 1].contourf(X_, Y_, U[:, 0, :])
# ax[2, 0].contourf(X_, Y_, u_[0, :, :])
# ax[2, 1].contourf(X_, Y_, U[0, :, :])
# fig.show()

# x, y = np.meshgrid(L, L)
# fig, ax = plt.subplots(2, 2)
# fig.colorbar(ax[0, 0].contourf(x, y, S[:, :, 0]), ax=ax[0, 0])
# fig.colorbar(ax[0, 1].contourf(x, y, s(x, y)), ax=ax[0, 1])
# fig.colorbar(ax[1, 0].contourf(x, y, u[:, :, 0]), ax=ax[1, 0])
# fig.colorbar(ax[1, 1].contourf(x, y, u(x, y)), ax=ax[1, 1])
# fig.show()

# plotting

fig, axes = plt.subplots(
    2, 2, sharex=True, sharey=False, figsize=(6.5, 4.5), dpi=100)

axes[0, 0].set_ylabel(r'Unfactored')
axes[1, 0].set_ylabel(r'Factored')
axes[0, 0].set_title(r'RMS Error')
axes[0, 1].set_title(r'Relative $\ell_\infty$ Error')

ax = axes[0, 0]
for Olim in Olims:
    ax.loglog(
        N3D, E[Olim]['rms'],
        label=common3d.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
# ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

ax = axes[0, 1]
for Olim in Olims:
    ax.loglog(
        N3D, E[Olim]['inf'],
        label=common3d.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
# ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

ax = axes[1, 0]
for Olim in Olims:
    ax.loglog(
        N3D, E_fac[Olim]['rms'],
        label=common3d.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
# ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

ax = axes[1, 1]
for Olim in Olims:
    ax.loglog(
        N3D, E_fac[Olim]['inf'],
        label=common3d.get_marcher_name(Olim).replace('_', ' '),
        linewidth=1, marker='|', markersize=3.5)
# ax.minorticks_off()
ax.set_xticks(N3D)
ax.set_xticklabels(map(str, N3D))

fig.legend()

fig.tight_layout()
fig.show()

# fig.savefig('../data/qv_plots_3d.eps')
