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

n = 10**3
x = np.linspace(-1., 1., n)
u = np.atan(40*x) + 2

def advect(row_size, width, advection_nodes):
    last_diff = 0.
    n_iter = 10**5
    advected = np.ones((n, row_size))
    for i in range(n_iter):
        diff = advected*u[:, np.newaxis]
        diff[ 1:, row_size//2:] = diff[:-1, row_size//2:] - diff[ 1:, row_size//2:]
        diff[:-1, :row_size//2] = diff[:-1, :row_size//2] - diff[ 1:, :row_size//2]
        diff[0, row_size//2:] *= 0
        diff[-1, :row_size//2] *= 0
        diff = 2e-1*(advection_nodes[np.newaxis, :]*diff + (1. - advected)/width*2/n)
        advected += diff
        last_diff = np.sqrt((diff*diff).sum())*n_iter
    print(f"Convergence error: {last_diff}")
    return advected

row_size = 6
width = .03

legendre_quad = legendre(row_size).weights
legendre_nodes = legendre_quad[:, 0]
legendre_weights = legendre_quad[:, 1]
proj_vec = legendre(row_size - 1)(legendre_nodes)*legendre_weights
"""
advected = advect(row_size, width, legendre_nodes)
for i_node in range(row_size):
    plt.plot(x, advected[:, i_node])
show()
proj = advected@proj_vec
plt.plot(x, proj)
show()
"""

row_size_advect = 4*row_size
width_smear = .5
x_smear = np.linspace(-1., 1., 100)
cheby_quad = chebyt(row_size_advect).weights
cheby_nodes = cheby_quad[:, 0]
advected = advect(row_size_advect, width, cheby_nodes)
fig, axs = plt.subplots(2, 1)
def compute_smoothness_mat(nodes):
    mat = np.zeros((len(nodes), row_size_advect))
    for i_node in range(len(nodes)):
        smoothness_nodes = legendre_nodes*width_smear*(nodes[i_node] + 1.)*(1. - nodes[i_node]) + nodes[i_node]
        interp_mat = interp(smoothness_nodes, cheby_nodes)
        mat[i_node, :] = proj_vec@interp_mat
    return mat
plot_smoothness = compute_smoothness_mat(x_smear)
smoothness_quad = legendre(3*row_size).weights
smoothness_nodes = smoothness_quad[:, 0]
smoothness_weights = smoothness_quad[:, 1]
quad_smoothness = compute_smoothness_mat(smoothness_nodes)
for plot_at in np.array([.1, .2, .4, .8])*width:
    plot_ind = int(n*plot_at) + n//2
    adv = advected[plot_ind, :]
    axs[0].scatter(cheby_nodes, adv, label=f"x = {plot_at:.3f}")
    axs[0].plot(x_smear, interp(x_smear, cheby_nodes)@adv)
    axs[1].plot(x_smear, plot_smoothness@adv)
    qs = quad_smoothness@adv
    axs[1].scatter(smoothness_nodes, qs, label=f"total = {np.sqrt((qs*smoothness_weights)@qs):.3e}")
axs[0].legend()
axs[1].legend()
for ax in axs:
    ax.grid(True)
fig.set_size_inches(20, 20)
fig, axs = plt.subplots(2, 1)
axs[0].plot(x, u)
smoothness_nodes = advected@quad_smoothness.transpose()
smoothness = np.sqrt(smoothness_nodes**2@smoothness_weights)
plt.plot(x, smoothness)
for ax in axs:
    ax.grid(True)
fig.set_size_inches(20, 20)
plt.show()
