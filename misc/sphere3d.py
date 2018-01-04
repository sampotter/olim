import sys

sys.path.insert(0, '../build/Release')

import eikonal as eik
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from skimage import measure

if __name__ == '__main__':
    n = 2**4 + 1
    m = eik.Olim26Rect(n, n, n, h=2/(n - 1), x0=1, y0=1, z0=1)
    i = int((n - 1)/2)
    m.addBoundaryNode(i, i, i)
    m.run()

    U = np.array(m)
    verts, faces, normals, values = measure.marching_cubes(U, 1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2])
    plt.show()
