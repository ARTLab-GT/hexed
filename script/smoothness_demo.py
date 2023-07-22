## \cond

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

def f(x):
    return np.tanh(1.*x)
x = np.linspace(-7, 7, int(1e3) + 1)
proj_width = 3

dx = x[1] - x[0]

def gradient(vec):
    grad = np.zeros(vec.size)
    grad[1:-1] = (vec[2:] - vec[:-2])/2/dx
    return grad

def diffusion(vec):
    res = np.zeros(vec.size)
    res[1:-1] = (vec[:-2] - 2*vec[1:-1] + vec[2:])/2/dx**2
    return res

#for n in [3]:
for n in range(1, 10):
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

    step = .9*dx**2
    sq = projd**2
    orig = sq + 0
    axs[0, 1].plot(x, orig)
    total = 0
    for j in range(3):
        diff = 4/rs
        pseudo = orig + 0
        iters = int(1e4)
        for it in range(iters):
            pseudo += step*(diffusion(pseudo) + (orig - pseudo)/diff)
            if it == iters//2:
                half = pseudo + 0
        orig = pseudo
        axs[0, 1].plot(x, half, linestyle="--")
        axs[0, 1].plot(x, pseudo)
        total += diff
    smeared = sq + 0
    for it in range(int(total/step)):
        smeared += step*diffusion(smeared)

    axs[0, 0].plot(x, shifted[n, :])
    axs[1, 0].plot(x, projd)
    axs[1, 1].plot(x, sq)
    axs[1, 1].plot(x, smeared)
    axs[1, 1].plot(x, half)

    for ax in axs.flatten():
        ax.grid(True)
    plt.show()

## \endcond
