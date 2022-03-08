from Basis import *
import numpy as np
import matplotlib.pyplot as plt
from sympy.integrals.quadrature import gauss_legendre, gauss_lobatto
from scipy.optimize import fsolve

calc_digits = 50
row_size = 6
nodes, weights = gauss_legendre(row_size, calc_digits)
nodes = [(node + 1)/2 for node in nodes]
weights = [weight/2 for weight in weights]
basis = Basis(nodes, weights, calc_digits=calc_digits)

n_elem = 8;
advection = np.zeros((n_elem*row_size, n_elem*row_size))
local_grad = np.zeros((n_elem*row_size, n_elem*row_size))
neighb_correction = np.zeros((n_elem*row_size, n_elem*row_size))
dissipation = np.zeros((n_elem*row_size, n_elem*row_size))
pos = np.zeros(n_elem*row_size)
global_weights = np.zeros(n_elem*row_size)
for i_elem in range(n_elem):
    for i_row in range(row_size):
        pos[i_elem*row_size + i_row] = i_elem + basis.nodes[i_row]
        global_weights[i_elem*row_size + i_row] = weights[i_row]
        for i_col in range(row_size):
            advection[  i_elem           *row_size + i_row, i_elem*row_size + i_col] = -basis.derivative(i_row, i_col)
            advection[  i_elem           *row_size + i_row, i_elem*row_size + i_col] -= basis.interpolate(i_col, 0)*basis.interpolate(i_row, 0)/weights[i_row]
            advection[((i_elem+1)%n_elem)*row_size + i_row, i_elem*row_size + i_col]  = basis.interpolate(i_col, 1)*basis.interpolate(i_row, 0)/weights[i_row]
            local_grad[ i_elem*row_size + i_row, i_elem*row_size + i_col] = basis.derivative(i_row, i_col)
            for i_side in [0, 1]:
                for j_side in [0, 1]:
                    boundary = basis.interpolate(i_col, 1-j_side)*basis.interpolate(i_row, 1-i_side)/weights[i_row]
                    neighb_correction[((i_elem+i_side)%n_elem)*row_size + i_row, ((i_elem+j_side)%n_elem)*row_size + i_col] += 0.5*(2*j_side - 1)*boundary
                    dissipation[      ((i_elem+i_side)%n_elem)*row_size + i_row, ((i_elem+j_side)%n_elem)*row_size + i_col] += (1 - 2*(j_side == i_side))*boundary
unstab_grad = local_grad + neighb_correction
ldg = unstab_grad@unstab_grad
ddg = unstab_grad@local_grad + dissipation
ident = np.identity(n_elem*row_size)
assert np.linalg.norm(global_weights@advection) < 1e-12
assert np.linalg.norm(global_weights@ldg) < 1e-12
assert np.linalg.norm(global_weights@ddg) < 1e-12

def euler(dt, diff_mat):
    return ident + dt*diff_mat
def rk(dt, step_weights, diff_mat):
    step_mat = ident
    for weight in step_weights:
        step_mat = (1. - weight)*ident + weight*euler(dt, diff_mat)@step_mat
    return step_mat
def ssp_rk3(dt, diff_mat):
    return rk(dt, [1., 1./4., 2./3.], diff_mat)

def max_cfl(matrix_fun, diff_mat):
    def error(log_dt):
        dt = np.exp(log_dt)
        eigvals, eigvecs = np.linalg.eig(matrix_fun(dt, diff_mat))
        return np.max(np.abs(eigvals)) - 1.
    log_dt = fsolve(error, 0.)
    dt = np.exp(log_dt)
    return dt[0], matrix_fun(dt, diff_mat)

mc, step_mat = max_cfl(ssp_rk3, advection)
eigvals, eigvecs = np.linalg.eig(step_mat)
plt.scatter(np.real(eigvals), np.imag(eigvals), label=f"Advection SSP RK3, CFL = {mc:.2g}")
mc, step_mat = max_cfl(euler, ldg)
eigvals, eigvecs = np.linalg.eig(step_mat)
plt.scatter(np.real(eigvals), np.imag(eigvals), label=f"LDG Fwd Euler, CFL = {mc:.2g}")
mc, step_mat = max_cfl(euler, ddg)
eigvals, eigvecs = np.linalg.eig(step_mat)
plt.scatter(np.real(eigvals), np.imag(eigvals), label=f"DDG Fwd Euler, CFL = {mc:.2g}")
angle = np.linspace(0, 2.*np.pi, 1000)
plt.plot(np.cos(angle), np.sin(angle), "k", label="stability boundary")
plt.title("eigenvalues")
plt.xlabel("real")
plt.xlim([plt.xlim()[0], 1.5])
plt.ylabel("imaginary")
plt.legend()
plt.gcf().set_size_inches(12, 8)
plt.show()
