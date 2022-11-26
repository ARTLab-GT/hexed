import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre, hermite, chebyt, jacobi

def interp(dest_nodes, src_nodes):
    mat = np.zeros((dest_nodes.size, src_nodes.size))
    for i_dest in range(dest_nodes.size):
        for i_src in range(src_nodes.size):
            lagrange = 1.
            for j_src in range(src_nodes.size):
                if i_src != j_src:
                    dest = dest_nodes[i_dest]
                    src0 = src_nodes[i_src]
                    src1 = src_nodes[j_src]
                    lagrange *= (dest - src1)/(src0 - src1)
            mat[i_dest, i_src] = lagrange
    return mat

def diffusion(vec):
    res = np.zeros(vec.size)
    res[1:-1] = vec[:-2] - 2*vec[1:-1] + vec[2:]
    return res

def f(x):
    return np.tanh(1.*x)
x = np.linspace(-7, 7, int(1e3) + 1)
proj_width = 3

#for n in range(1, 10):
for n in [3]:
    rs = 2*n
    quad = legendre(rs + 1).weights
    nodes = quad[:, 0]
    weights = quad[:, 1]
    fig, axs = plt.subplots(2, 2)
    fig.set_size_inches(14, 8)

    poly = legendre(rs)
    shifted = f(proj_width*nodes[:, np.newaxis] + x)
    proj = poly(nodes)*weights
    projd = proj@shifted
    smeared = projd**2
    for it in range(10000):
        smeared += .45*diffusion(smeared)

    axs[0, 0].plot(x, shifted[n, :])
    axs[1, 0].plot(x, projd)
    axs[1, 1].plot(x, projd**2)
    axs[1, 1].plot(x, smeared)

    for ax in axs.flatten():
        ax.grid(True)
    plt.show()
