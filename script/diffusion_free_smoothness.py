import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre, hermite, chebyt, jacobi

row_size = 6
width = .03

alpha0 = 1.
alpha1 = 0.
plt.scatter(jacobi(row_size, alpha0, alpha1).weights[:, 0], np.zeros(row_size))
plt.scatter(jacobi(row_size, alpha1, alpha0).weights[:, 0], np.zeros(row_size))
plt.grid(True)
plt.show()

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

legendre_quad = legendre(row_size).weights
legendre_nodes = legendre_quad[:, 0]
legendre_weights = legendre_quad[:, 1]

advected = advect(row_size, width, legendre_nodes)
for degree in range(row_size):
    proj_vec = width**((row_size - degree)/2)*legendre(degree)(legendre_nodes)*legendre_weights
    proj = advected@proj_vec
    plt.plot(x, proj)
plt.grid(True)
plt.gcf().set_size_inches(20, 10)
plt.show()
exit()

row_size_advect = 3*row_size
x_smear = np.linspace(-1., 1., 100)
cheby_quad = chebyt(row_size_advect).weights
cheby_nodes = cheby_quad[:, 0]
advected = advect(row_size_advect, width, cheby_nodes)
fig, axs = plt.subplots(2, 1)
def compute_smoothness_mat(coefs):
    mat = np.zeros((coefs.shape[0], row_size_advect))
    for i in range(coefs.shape[0]):
        jacobi_quad = jacobi(row_size, coefs[i, 0], coefs[i, 1]).weights
        interp_mat = interp(jacobi_quad[:, 0], cheby_nodes)
        mat[i, :] = (jacobi(row_size - 1, coefs[i, 0], coefs[i, 1])(jacobi_quad[:, 0])*jacobi_quad[:, 1])@interp_mat
    return mat
alpha0 = 0.
alpha1 = 6.
plot_coefs = np.array([np.linspace(alpha0, alpha1, 100), np.linspace(alpha1, alpha0, 100)]).transpose()
plot_smoothness = compute_smoothness_mat(plot_coefs)
transformed_nodes = .5*(legendre_nodes + 1)
quad_coefs = np.array([alpha0 + (alpha1 - alpha0)*transformed_nodes,
                       alpha1 + (alpha0 - alpha1)*transformed_nodes]).transpose()
quad_smoothness = compute_smoothness_mat(quad_coefs)
smoothness_weights = legendre_weights*(1 + legendre_nodes)*(1 - legendre_nodes)
for plot_at in np.array([.1, .2, .4, .8, 1.6])*width:
    plot_ind = int(n*plot_at) + n//2
    adv = advected[plot_ind, :]
    axs[0].scatter(cheby_nodes, adv, label=f"x = {plot_at:.3f}")
    axs[0].plot(x_smear, interp(x_smear, cheby_nodes)@adv)
    axs[1].plot(plot_coefs[:, 0], plot_smoothness@adv)
    axs[1].scatter(quad_coefs[:, 0], quad_smoothness@adv)
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
