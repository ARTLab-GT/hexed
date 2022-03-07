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

n_elem = 4;
diff_mat = np.zeros((n_elem*row_size*2, n_elem*row_size*2)) # time derivative matrix
pos = np.zeros(n_elem*row_size)
global_weights = np.zeros((2, 2*n_elem*row_size))
flux_jacobian = np.array([[1, 2], [3, 4]])
char_speeds = np.linalg.eig(flux_jacobian)[0]
mcs = np.max(char_speeds)
for i_var in range(2):
    for j_var in range(2):
        var_mat = diff_mat[j_var*n_elem*row_size:(j_var+1)*n_elem*row_size, j_var*n_elem*row_size:(j_var+1)*n_elem*row_size]
        for i_elem in range(n_elem):
            for i_row in range(row_size):
                pos[i_elem*row_size + i_row] = i_elem + basis.nodes[i_row]
                global_weights[i_var, i_var*n_elem*row_size + i_elem*row_size + i_row] = weights[i_row]
                for i_col in range(row_size):
                    var_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] = -basis.derivative(i_row, i_col)*flux_jacobian[i_var, j_var]
                    var_mat[  i_elem           *row_size + i_row, i_elem*row_size + i_col] -= basis.interpolate(i_col, 0)*basis.interpolate(i_row, 0)/weights[i_row]
                    var_mat[((i_elem+1)%n_elem)*row_size + i_row, i_elem*row_size + i_col]  = basis.interpolate(i_col, 1)*basis.interpolate(i_row, 0)/weights[i_row]
plt.imshow(diff_mat)
plt.show()
ident = np.identity(2*n_elem*row_size)
assert np.linalg.norm(global_weights@diff_mat) < 1e-12

def euler(dt):
    return ident + dt*diff_mat/mcs
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
        return np.min(np.real(eigvals)) + 1.
    log_dt = fsolve(error, 0.)
    return np.exp(log_dt)

mc = max_cfl(euler)
eigvals, eigvecs = np.linalg.eig(euler(mc))
plt.scatter(np.real(eigvals), np.imag(eigvals))
plt.show()
print(max_cfl(euler))
print(max_cfl(ssp_rk3))
