from Basis import *
import numpy as np
import matplotlib.pyplot as plt
from sympy.integrals.quadrature import gauss_legendre, gauss_lobatto
from scipy.optimize import fsolve

calc_digits = 50
row_size = 6
nodes, weights = gauss_legendre(row_size, calc_digits)
#nodes, weights = gauss_lobatto(row_size, calc_digits)
nodes = [(node + 1)/2 for node in nodes]
weights = [weight/2 for weight in weights]
basis = Basis(nodes, weights, calc_digits=calc_digits)

n_elem = 2;
diff_mat = np.zeros((n_elem*row_size, n_elem*row_size)) # time derivative matrix
pos = np.zeros(n_elem*row_size)
global_weights = np.zeros(n_elem*row_size)
for i_elem in range(n_elem):
    for i_row in range(row_size):
        pos[i_elem*row_size + i_row] = i_elem + basis.nodes[i_row]
        global_weights[i_elem*row_size + i_row] = weights[i_row]
        for i_col in range(row_size):
            diff_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] = -basis.derivative(i_row, i_col)
            diff_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] -= basis.interpolate(i_col, 0)*basis.interpolate(i_row, 0)/weights[i_row]
            diff_mat[((i_elem+1)%n_elem)*row_size + i_row, i_elem*row_size + i_col]  = basis.interpolate(i_col, 1)*basis.interpolate(i_row, 0)/weights[i_row]
ident = np.identity(n_elem*row_size)
assert np.linalg.norm(global_weights@diff_mat) < 1e-12

def euler(dt):
    return ident + dt*diff_mat
def rk(dt, step_weights):
    step_mat = ident
    for weight in step_weights:
        step_mat = (1. - weight)*ident + weight*euler(dt)@step_mat
    return step_mat
def ssp_rk3(dt):
    return rk(dt, [1., 1./4., 2./3.])

def max_cfl(matrix_fun):
    def error(log_dt):
        dt = np.exp(log_dt)
        eigvals, eigvecs = np.linalg.eig(matrix_fun(dt))
        return np.max(np.abs(eigvals)) - 1.
    log_dt = fsolve(error, 0.)
    return np.exp(log_dt)

mc = max_cfl(ssp_rk3)[0]
print(mc)
eigvals, eigvecs = np.linalg.eig(ssp_rk3(mc))
plt.scatter(np.real(eigvals), np.imag(eigvals))
angle = np.linspace(0, 2.*np.pi, 1000)
plt.plot(np.cos(angle), np.sin(angle), "k")
plt.show()
