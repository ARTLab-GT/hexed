import sympy as sp
import numpy as np
from scipy.optimize import fsolve, minimize
import warnings

## \namespace Basis \brief module for `Basis.Basis`

class Basis:
    r"""!
    Computes numerical parameters for nodal polynmial bases (such as hexed::Gauss_legendre and hexed::Gauss_lobatto)
    that are based on [Gaussian quadrature rules](https://en.wikipedia.org/wiki/Gaussian_quadrature).
    Used in `auto_generate.py`.
    Calculations for most functions are performed in arbitrary-precision arithmetic with sympy
    to minimize roundoff errors (which can sometimes be nonnegligible for high-order bases
    if the calculations are performed with the same precision as the output).
    """

    def __init__(self, nodes, weights, repr_digits=20, calc_digits=50):
        r"""!
        \param nodes Quadrature nodes. Should be in interval [0, 1]
        \param weights Quadrature weights. Should sum to 1.
        \param repr_digits number of digits returned in output
        \param calc_digits number of digits used in calculations (recommended to be greater than `repr_digits`
        \attention Standard convention for Gaussian quadratures
        is that nodes are in [-1, 1] and weights sum to 2.
        Be sure to perform any necessary conversions to ensure that nodes are in [0, 1] and weights sum to 1.
        """
        ## \private
        self.repr_digits = repr_digits
        ## \private
        self.calc_digits = calc_digits
        ## \private
        self.row_size = len(nodes)
        assert len(weights) == self.row_size
        ## \private
        self.nodes = [sp.Float(node, self.calc_digits) for node in nodes]
        ## \private
        self.weights = [sp.Float(weight, self.calc_digits) for weight in weights]
        ## \private
        self.ortho = []
        for i in range(self.row_size):
            self.ortho.append(self.legendre(i))

    def node(self, i):
        r"""! \brief returns the `i`th quadrature node
        \note arbitrary precision
        """
        return sp.Float(self.nodes[i], self.repr_digits)

    def weight(self, i):
        r"""! \brief returns the `i`th quadrature weight
        \note arbitrary precision
        """
        return sp.Float(self.weights[i], self.repr_digits)

    def derivative(self, i_result, i_operand):
        r"""! \brief geta an element of the differentiation matrix.
        \details Computes the derivative of the `i_operand`th basis polynomial at the `i_result`th node
        \note arbitrary precision
        """
        if i_result == i_operand:
            dp = sp.Float(0, self.calc_digits)
            for n in range(self.row_size):
                if n != i_result:
                    dp += sp.Float(1, self.calc_digits)/(self.nodes[i_result] - self.nodes[n])
        else:
            dp = sp.Float(1, self.calc_digits)
            for n in range(self.row_size):
                if n != i_operand:
                    dp /= sp.Float(self.nodes[i_operand] - self.nodes[n])
                    if n != i_result:
                        dp *= sp.Float(self.nodes[i_result] - self.nodes[n])
        return sp.Float(dp, self.repr_digits)

    def interpolate(self, i, position, calc=False):
        r"""! \brief interpolates one of the basis polynomials to a specified point
        \param i will interpolate the polynomial associated with the `i`th node
        \param position position to interpolate to
        \param calc if True, will return `self.calc_digits` digits. Otherwise, will return `self.repr_digits` digits.
        \note arbitrary precision
        """
        pos = sp.Float(position, self.calc_digits)
        nodes = self.nodes.copy()
        main_node = nodes.pop(i)
        result = sp.Float(1, self.calc_digits)
        for node in nodes:
            result *= (pos - node)/(main_node - node)
        return sp.Float(result, self.calc_digits if calc else self.repr_digits)

    def legendre(self, degree):
        r"""! \brief compute a Legendre polynomial
        \param degree degree of the polynomial to compute
        \returns a list of the values of the Legendre polynomial at each node
        \details Polynomial is has norm 1 with respect to the quadrature under consideration.
        That is, the sum of the squares of the returned values multiplied by the quadrature weights is 1.
        \note arbitrary precision
        """
        x = sp.Symbol("x")
        poly = sp.legendre(degree, x)
        vals = []
        norm = sp.Float(0, self.calc_digits)
        for i_node in range(self.row_size):
            node = self.nodes[i_node]
            vals.append(poly.subs(x, node*2 - 1).evalf(self.calc_digits))
            norm += vals[i_node]**2*self.weights[i_node]
        norm = norm**sp.Rational(1, 2)
        for i_node in range(self.row_size):
            vals[i_node] /= norm
        return vals

    def get_ortho(self, degree, i_node):
        r"""! equivalent to `self.legendre(degree)[i_node]` """
        return sp.Float(self.ortho[degree][i_node], self.repr_digits)

    def prolong(self, i_result, i_operand, i_half, calc=False):
        r"""! \brief gets an element of the prolongation matrix
        \details Suppose you split the interval [0, 1] in half
        and defined another Basis on each half with the same quadrature rule.
        Then this function evaluates the `i_result`th basis polynomial at the `i_operand`th quadrature point
        of the `i_half`th fine (half-interval) basis.
        Parameter `calc` is forwarded to `interpolate`.
        \note arbitrary precision
        """
        position = (self.nodes[i_result] + i_half)/2
        return self.interpolate(i_operand, position, calc)

    def restrict(self, i_result, i_operand, i_half):
        r"""! \brief get an element of the restriction matrix
        \details pseudo-inverse of prolongation matrix computed by `prolong`.
        \f$L_2\f$ projection of `i_operand`th polynomial on the `i_half`th fine basis (see `prolong`)
        with `i_result`th  polynomial of this basis, computed by the quadrature rule.
        \attention Doesn't work for Gauss-Lobatto quadrature
        \note arbitrary precision
        """
        result = self.weights[i_operand]/2*self.prolong(i_operand, i_result, i_half, True)/self.weights[i_result]
        return sp.Float(result, self.repr_digits)

    def time_coefs(self, n_elem = 16):
        r"""! \brief Compute coefficients for optimized 2-stage time integration scheme.
        \param n_elem number of elements to use
        \details Based on numerical eigenvalue calculations for the linear advection and diffusion equations
        on a 1D mesh with periodic boundary conditions.
        \returns a pair of pairs:
        - `time_coefs()[0][0]` is the maximum stable CFL number for advection
        - `time_coefs()[0][1]` is the cancellation coefficient for advection
        - `time_coefs()[1][0]` is the maximum stable CFL number for diffusion
        - `time_coefs()[1][1]` is the cancellation coefficient for diffusion
        \note standard floating-point precision (whatever that is for Python -- I think double)
        """
        # sorry for the lack of comments... remind me to get back to this later
        global_weights = np.zeros(n_elem*self.row_size)
        local_grad = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_avrg = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_jump = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        global_nodes = global_weights*0
        for i_elem in range(n_elem):
            for i_row in range(self.row_size):
                global_weights[i_elem*self.row_size + i_row] = self.weights[i_row]
                global_nodes[i_elem*self.row_size + i_row] = i_elem + self.nodes[i_row]
                for i_col in range(self.row_size):
                    row = i_elem*self.row_size + i_row
                    col = i_elem*self.row_size + i_col
                    local_grad[row, col] = self.derivative(i_row, i_col)
                    for i_side in [0, 1]:
                        for j_side in [0, 1]:
                            row = ((i_elem+i_side)%n_elem)*self.row_size + i_row
                            col = ((i_elem+j_side)%n_elem)*self.row_size + i_col
                            boundary = self.interpolate(i_col, 1-j_side)*self.interpolate(i_row, 1-i_side)/self.weights[i_row]
                            neighb_avrg[row, col] += 0.5*(2*j_side - 1)*boundary
                            neighb_jump[row, col] += 0.5*(1 - 2*(j_side == i_side))*boundary
        advection = -local_grad + -neighb_avrg + neighb_jump
        diffusion = (local_grad + neighb_avrg)@(local_grad + neighb_avrg)
        proj = np.zeros((self.row_size, self.row_size))
        interp = np.zeros((self.row_size, self.row_size))
        for row in range(self.row_size):
            for col in range(self.row_size):
                proj[row, col] = self.get_ortho(row, col)*self.weights[col]
                interp[row, col] = self.get_ortho(col, row)
        assert np.linalg.norm(interp@proj - np.identity(self.row_size)) < 1e-12
        scale = .5**np.arange(self.row_size)
        interp *= scale
        filt = interp@proj
        ident = np.identity(n_elem*self.row_size)
        assert np.linalg.norm(global_weights@advection) < 1e-12
        assert np.linalg.norm(global_weights@diffusion) < 1e-12
        """
        filtered = advection + 0
        for i_elem in range(n_elem):
            filtered[i_elem*self.row_size:(i_elem+1)*self.row_size, :] = filt@filtered[i_elem*self.row_size:(i_elem+1)*self.row_size, :]
        safety = 0.8
        max_cfl = np.real(-2*safety/np.linalg.eig(advection)[0].min())
        print(max_cfl)
        step = max_cfl*advection
        eigvals, eigvecs = np.linalg.eig(ident + step)
        plt.scatter(np.real(eigvals), np.imag(eigvals), marker="+")
        eigvals, eigvecs = np.linalg.eig(ident + step + .5/safety*step@step)
        plt.scatter(np.real(eigvals), np.imag(eigvals), marker="+")
        max_cfl = np.real(-2*safety/np.linalg.eig(filtered)[0].min())
        print(max_cfl)
        step = max_cfl*filtered
        eigvals, eigvecs = np.linalg.eig(ident + step)
        plt.scatter(np.real(eigvals), np.imag(eigvals), marker="+")
        eigvals, eigvecs = np.linalg.eig(ident + step + .5/safety*step@step)
        plt.scatter(np.real(eigvals), np.imag(eigvals), marker="+")
        angle = np.linspace(-np.pi, np.pi, 1000)
        plt.plot(np.cos(angle), np.sin(angle), color="k")
        angle = np.linspace(-1, 1, 1000)
        plt.plot(1 - 0.5/safety*angle**2, angle, color="k", linestyle="dashed")
        plt.axis("equal")
        plt.grid(True)
        plt.show()
        assert False
        """
        class polynomial:
            def __init__(self, coefs):
                ## \private
                self.coefs = coefs
            def __call__(self, dt, mat):
                step_mat = ident + dt*mat
                mat_pow = mat + 0
                for coef in self.coefs:
                    mat_pow = mat @ mat_pow
                    step_mat += dt*coef*mat_pow
                return step_mat
        def max_cfl(mat, time_scheme):
            "numerically estimates the maximum CFL number of a scheme"
            def error(log_dt):
                dt = np.exp(log_dt)
                eigvals, eigvecs = np.linalg.eig(time_scheme(dt, mat))
                return np.max(np.abs(eigvals)) - 1.
            return np.exp(fsolve(error, -np.log(self.row_size))[0])
        if True:
            def euler(dt, mat):
                return ident + dt*mat
            def rk(dt, step_weights, mat):
                step_mat = ident
                for weight in step_weights:
                    step_mat = (1. - weight)*ident + weight*euler(dt, mat)@step_mat
                return step_mat
            def ssp_rk3(dt, mat):
                return rk(dt, [1., 1./4., 2./3.], mat)
        safety = 0.9
        cfl_adv = -2*safety/np.linalg.eig(advection)[0].real.min()
        min_eig = np.linalg.eig(diffusion)[0].real.min()
        cfl_diff = -8*safety/min_eig
        coefs = ((cfl_adv, .5/safety*cfl_adv), (cfl_diff, -1/min_eig))
        print(f"        {self.row_size - 1} & {coefs[0][0]:.4f} & {coefs[0][1]:.5f} & {coefs[1][0]:.4f} & {coefs[1][1]:.5f} & {max_cfl(advection, ssp_rk3):.5f} & {max_cfl(diffusion, ssp_rk3):.5f} \\\\")
        return coefs

"""
from sympy.integrals.quadrature import gauss_legendre
import matplotlib.pyplot as plt
nodes, weights = gauss_legendre(6, 50)
nodes = [(node + 1)/2 for node in nodes]
weights = [weight/2 for weight in weights]
basis = Basis(nodes, weights, calc_digits=50)
basis.time_coefs()
"""
