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

#for n in range(1, 10):
for n in [3]:
    rs = 2*n
    quad = legendre(rs + 1).weights
    nodes = quad[:, 0]
    weights = quad[:, 1]
    fig, axs = plt.subplots(2, 2)
    fig.set_size_inches(14, 8)

    smear_rs = rs
    #scale = 10/rs
    scale = proj_width
    """
    smear_weight = np.exp(-(x/scale)**2)
    smear_quad = hermite(2*smear_rs + 1).weights
    """
    """
    smear_weight = x*0
    smear_weight[np.logical_and(x > -scale, x < scale)] = 1
    smear_quad = legendre(2*smear_rs + 1).weights
    """
    power = 3
    smear_weight = x*0
    mask = np.logical_and(x > -scale, x < scale)
    smear_weight[mask] = ((scale - x[mask])*(scale + x[mask])/scale/smear_rs)**power
    smear_quad = jacobi(2*smear_rs - 1, power, power).weights
    smear_quad[:, 0] *= scale
    """
    shift_rs = smear_rs + 3
    smear_shift = chebyt(shift_rs).weights[:, 0]*scale
    """
    shift_rs = smear_quad.shape[0]
    smear_shift = smear_quad[:, 0]

    poly = legendre(rs)
    """
    sample_width = scale + proj_width
    shift = chebyt(30).weights[:, 0]*sample_width
    raw = f(shift[:, np.newaxis] + x)
    shifted = np.zeros((smear_quad.shape[0], rs + 1, x.size))
    for i_node in range(rs + 1):
        shifted[:, i_node, :] = interp(smear_quad[:, 0] + proj_width*nodes[i_node], shift)@raw
    """
    shifted = f(proj_width*nodes[np.newaxis, :, np.newaxis] + smear_shift[:, np.newaxis, np.newaxis] + x)
    proj = poly(nodes)*weights
    projd = proj@shifted
    smeared = smear_quad[:, 1]@interp(smear_quad[:, 0], smear_shift)@projd**2

    axs[0, 0].plot(x, shifted[shift_rs//2, n, :])
    axs[1, 0].plot(x, projd[shift_rs//2])
    axs[0, 1].plot(x, smear_weight)
    axs[0, 1].scatter(smear_quad[:, 0], smear_quad[:, 1])
    axs[1, 1].plot(x, projd[shift_rs//2]**2)
    axs[1, 1].plot(x, smeared)

    for ax in axs.flatten():
        ax.grid(True)
    plt.show()
