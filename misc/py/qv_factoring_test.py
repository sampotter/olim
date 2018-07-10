import sys; sys.path.insert(0, '../../build/Release')
import matplotlib.pyplot as plt; plt.ion()
import numpy as np
import pyeikonal as eik

from matplotlib.colors import LogNorm

vx, vy = 5, 20

s = lambda x, y: 1/(2 + vx*x + vy*y)

def make_u(x_fac, y_fac, vx, vy, s_fac, s):
    return lambda x, y: \
        (1/np.sqrt(vx**2 + vy**2)) * \
        np.arccosh(
            1 +
            s_fac*s(x, y)*(vx**2 + vy**2)*((x - x_fac)**2 + (y - y_fac)**2)/2)

x_fac_1, y_fac_1 = 0, 0
x_fac_2, y_fac_2 = 0.8, 0

s_fac_1 = s(x_fac_1, y_fac_1)
s_fac_2 = s(x_fac_2, y_fac_2)

u_1 = make_u(x_fac_1, y_fac_1, vx, vy, s_fac_1, s)
u_2 = make_u(x_fac_2, y_fac_2, vx, vy, s_fac_2, s)

u = lambda x, y: np.minimum(u_1(x, y), u_2(x, y))

Olim = eik.Olim8Rect

R_fac = 2

N = 2**np.arange(5, 14) + 1

E = {'2': np.zeros(len(N)), 'inf': np.zeros(len(N))}
E_fac = {'2': np.zeros(len(N)), 'inf': np.zeros(len(N))}

for k, n in enumerate(N):
    print('n = %d (%d/%d)' % (n, k + 1, len(N)))

    L = np.linspace(0, 1, n)
    X, Y = np.meshgrid(L, L)
    S = s(X, Y)

    h = 1/(n - 1)
    i_1, j_1 = int(y_fac_1/h), int(x_fac_1/h)
    i_2, j_2 = int(y_fac_2/h), int(x_fac_2/h)
    R_fac = 0.1

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
    U_fac = np.array([[m_fac.getValue(i, j) for j in range(n)] for i in range(n)])

    E['2'][k] = np.linalg.norm(U - u(X, Y))/np.linalg.norm(U)
    E['inf'][k] = np.linalg.norm(U - u(X, Y), np.inf)/np.linalg.norm(U, np.inf)

    E_fac['2'][k] = np.linalg.norm(U_fac - u(X, Y))/np.linalg.norm(U_fac)
    E_fac['inf'][k] = \
        np.linalg.norm(U_fac - u(X, Y), np.inf)/np.linalg.norm(U_fac, np.inf)

fig = plt.figure()

ax = fig.add_subplot(3, 2, 1)
ax.imshow(u(X, Y))
ax.invert_yaxis()

ax = fig.add_subplot(3, 2, 2)
ax.imshow(U)
ax.invert_yaxis()

eps = np.finfo(np.float64).eps
err = np.abs(u(X, Y) - U)
ax = fig.add_subplot(3, 2, 3)
ax.imshow(np.maximum(eps, err), norm=LogNorm(vmin=eps, vmax=err.max()))
ax.invert_yaxis()

eps = np.finfo(np.float64).eps
err = np.abs(u(X, Y) - U_fac)
ax = fig.add_subplot(3, 2, 4)
ax.imshow(np.maximum(eps, err), norm=LogNorm(vmin=eps, vmax=err.max()))
ax.invert_yaxis()

ax = fig.add_subplot(3, 2, 5)
ax.loglog(N, E['2'], label='E')
ax.loglog(N, E_fac['2'], label='E_fac')
ax.legend()

ax = fig.add_subplot(3, 2, 6)
ax.loglog(N, E['inf'], label='E')
ax.loglog(N, E_fac['inf'], label='E_fac')
ax.legend()

fig.show()

R_1 = np.sqrt((X - x_fac_1)**2 + (Y - y_fac_1)**2)
R_2 = np.sqrt((X - x_fac_2)**2 + (Y - y_fac_2)**2)

I_1, J_1 = np.where(R_1 <= R_fac)
I_2, J_2 = np.where(R_2 <= R_fac)

