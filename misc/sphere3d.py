import sys

sys.path.insert(0, '../build/Release')

import eikonal as eik
import matplotlib.pyplot as plt
import numpy as np
import skimage

if __name__ == '__main__':
    n = 2**5 + 1
    m = eik.BasicMarcher3D(n, n, n, h=2/(n - 1), x0=1, y0=1, z0=1)
    i = int((n - 1)/2)
    m.addBoundaryNode(i, i, i)
    m.run()
    U = np.array(m)
    verts, faces = skimage.measure.marching_cubes(U, 1)
    
