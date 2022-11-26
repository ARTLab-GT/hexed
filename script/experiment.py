import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre, hermite, chebyt

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
plt.plot(x, np.exp(x/7), color="k")
nodes = chebyt(4).weights[:, 0]*7
plt.scatter(nodes, np.exp(nodes/7))
plt.plot(x, interp(x, nodes)@np.exp(nodes/7))
plt.show()
exit()

for n in range(1, 10):
#for n in [3]:
    rs = 2*n

    quad = legendre(rs + 1).weights
    nodes = quad[:, 0]
    weights = quad[:, 1]
    fig, axs = plt.subplots(2, 2)
    fig.set_size_inches(14, 8)

    smear_rs = 50
    scale = 10/rs
    smear_weight = np.exp(-(x/scale)**2)
    smear_quad = hermite(2*smear_rs + 1).weights
    smear_quad[:, 0] *= scale

    poly = legendre(rs)
    shifted = f(3*nodes[np.newaxis, :, np.newaxis] + smear_quad[:, 0, np.newaxis, np.newaxis] + x)
    proj = poly(nodes)*weights
    projd = proj@shifted
    smeared = smear_quad[:, 1]@projd**2

    axs[0, 0].plot(x, shifted[smear_rs, n, :])
    axs[1, 0].plot(x, projd[smear_rs])
    axs[0, 1].plot(x, smear_weight)
    axs[0, 1].scatter(smear_quad[:, 0], smear_quad[:, 1])
    axs[1, 1].plot(x, projd[smear_rs]**2)
    axs[1, 1].plot(x, smeared)

    for ax in axs.flatten():
        ax.grid(True)
    plt.show()
