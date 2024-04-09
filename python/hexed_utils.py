import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import time
import re
import pandas as pd
import sympy as sp
from scipy.optimize import fsolve, minimize
import warnings
from sympy.integrals.quadrature import gauss_legendre, gauss_lobatto

## \namespace hexed_utils
# \brief A Python module with random useful tools.

def naca(desig, n_points = 1000, closure = "warp"):
    r"""! \brief Constructs a NACA 4-digit airfoil geometry.
    \details Returns an n by 2 numpy array representing the coordinates of the airfoil at discrete points.
    This array can then be passed to `Solver.generate_mesh` as a geometry.
    Points are clustered near the leading edge but not the trailing (see implementation for details).
    This function is the recommended way to generate NACA airfoil geometry for Hexed simulations,
    as importing airfoils from coordinate files requires some \ref geom_fitting "special care".
    \param desig String representing the airfoil designation (e.g., `"0012"` for the NACA0012).
                 In general, we cannot accept this parameter as an `int` because of possible leading zeros.
    \param n_points Number of points on the airfoil surface. Don't be stingy, since DG is finnicky with discrete geometry representations --
                    1000 is actually on the lower end of what I normally use.
    \param closure If and how to close the trailing edge. There are 3 options:
                   - `"warp"`: Close the trailing edge by adding a 4th-degree polyomial of \f$ x_0 \f$ to \f$ x_1 \f$.
                   - `"segment"`: Close the trailing edge by adding a line segment connecting the last point to the first point,
                     causing the array to be `(n_points + 1)*2` instead of `n_points*2`.
                   - `"none"`, `None`, or `False`: Don't close the trailing edge.
    """
    try:
        desig = str(desig)
        assert len(desig) == 4
    except Exception as e:
        raise User_error("cannot interpret `desig` as a 4-character string") from e
    camber_max = int(desig[0])
    camber_loc = int(desig[1])
    thickness = int(desig[2:])*1e-2
    coords = np.zeros((n_points, 2))
    param = np.linspace(-1., 1., n_points)
    coords[:, 0] = param**2
    ap = np.abs(param)
    coords[:, 1] = 5*thickness*param*(.2969 - .1260*ap - .3516*ap**3 + .2843*ap**5 - .1015*ap**7)
    if camber_loc > 0:
        camber_max *= 1e-2
        camber_loc *= 1e-1
        s = coords[:, 0] <  camber_loc
        coords[s, 1] += camber_max/camber_loc**2*(2*camber_loc*param[s]**2 - param[s]**4)
        s = coords[:, 0] >= camber_loc
        coords[s, 1] += camber_max/(1 - camber_loc)**2*(1 - 2*camber_loc + 2*camber_loc*param[s]**2 - param[s]**4)
    if closure == "warp":
        coords[:, 1] -= param*ap**7*coords[-1, 1]
    elif closure == "segment":
        coords = np.concatenate([coords, coords[[0], :]])
    elif closure and closure.lower() != "none":
        raise User_error("unrecognized value of `closure` parameter")
    return coords

def joukowsky(thickness, camber = 0., n_points = 1000, scale = True):
    r"""! \brief constructs a [Joukowsky airfoil](https://en.wikipedia.org/wiki/Joukowsky_transform)
    \details This is a family of airfoils with a cusped trailing edge which have analytic solutions for the incompressible flow around them.
    \param thickness _approximate_ thickness-to-chord ratio
    \param camber (radian) angle between the trailing edge and the chord line
    \param n_points Number of points on the airfoil surface
    \param scale If `True`, scale the airfoil so that the leading edge is at (0, 0) and the trailing edge is at (0, 1).
        Otherwise, the airfoil is left at the size and position dictated by the Joukowsky transform.
    """
    # create circle in complex plain
    points = np.exp(np.linspace(0, 2*np.pi, n_points)*1j)
    # transform circle so that it has the correct position relative to the singularity of the Joukowsky transform
    radius = 1 + 4*thickness/3**1.5
    points *= radius
    points = (points - radius)*np.exp(camber*-.5j) + 1
    # apply Joukowsky transform
    points = points + 1/points
    # scale airfoil to match engineering conventions
    if scale:
        points -= points.real.min()
        points /= points.real.max()
    return np.array([points.real, points.imag]).transpose()

## \brief alternative transliteration
zhukovsky = joukowsky

class History_plot:
    r"""! \brief creates a real-time, interactive plot of the convergence history
    \details Convergence history is obtained from the `output.txt` file which contains the console output of \ref hexecute.
    Every column in the output whose name is not in `History_plot.column_blacklist` will be plotted in its own subplot.
    Columns with names ending in `residual` will be plotted on a log scale.
    A `History_plot` instance should be created in a separate process from the solver
    but may be created before, after, or during the simulation---the plot will not appear until an `output.txt` file exists.
    A `History_plot` can be created directly from the \ref hil "HIL" solver script with $\ref plot_history.
    """
    ## \brief names of output columns __not__ to plot
    ## \details You may modify this variable for the class or for instances
    column_blacklist = ["flow_time", "time_step"]

    def _infinite_generator(self):
        while not self._stop:
            yield None

    def _read_new_lines(self):
        if not os.path.exists(self._directory + "output.txt"): return
        with open(self._directory + "output.txt", "rb") as output_file:
            output_file.seek(self._file_position, 0)
            self._lines += [line.decode("utf-8") for line in output_file.readlines()]
            self._file_position = output_file.tell()

    def __init__(self, directory = "hexed_out", interval = 0.2):
        r"""! \brief creates and shows an animated history plot
        \param directory Convergence history will be obtained by looking for a file `output.txt` in `directory`.
        \param interval The plot will be updated every `interval` seconds to include new iterations.
        """
        if directory[-1] != "/": directory += "/"
        self._directory = directory
        self._interval = interval
        self._file_position = 0
        self._data = None
        self._lines = []
        while self._data is None:
            self._read_new_lines()
            while self._lines:
                line = self._lines.pop(0).replace(" ", "").replace("\n", "")
                if re.match("iteration,", line):
                    self._data = pd.DataFrame(columns = line.split(","))
                    break
            time.sleep(self._interval)
        self._plot_columns = [col for col in self._data.columns if col not in self.column_blacklist + ["iteration"]]
        self._stop = False
        self._fig, self._axs = plt.subplots(1, len(self._plot_columns))
        plt.tight_layout()
        self._fig.set_size_inches(18, 5)
        ani = FuncAnimation(self._fig, self._update, frames = self._infinite_generator, init_func = self._init,
                            blit = True, repeat = False, interval = int(self._interval*1e3), cache_frame_data = False)
        plt.show()

    def _init(self):
        self._curves = []
        for i_col in range(len(self._plot_columns)):
            col = self._plot_columns[i_col]
            ax = self._axs[i_col]
            self._curves.append(ax.plot([], [])[0])
            ax.set_xlim(0., 1.)
            ax.grid(True)
            ax.set_xlabel("iteration")
            ax.set_ylabel(col)
            if col.endswith("residual"):
                self._axs[i_col].set_ylim(0.1, 1.)
                self._axs[i_col].set_yscale("log")
        return self._curves

    def _update(self, _):
        self._read_new_lines()
        while self._lines:
            line = self._lines.pop(0)
            if line.startswith("simulation complete"):
                self._stop = True
            if re.match(" *[0-9]+,", line):
                entries = line.split(",")
                self._data.loc[self._data.shape[0]] = [int(entries[0])] + [float(e) for e in entries[1:]]
                last_iter = self._data["iteration"][self._data.shape[0] - 1]
                if last_iter > self._axs[0].get_xlim()[1]:
                    for ax in self._axs:
                        ax.set_xlim(0, ax.get_xlim()[1]*2)
            for i_col in range(len(self._plot_columns)):
                ax = self._axs[i_col]
                col = self._plot_columns[i_col]
                last_value = self._data[col][self._data.shape[0] - 1]
                if col.endswith("residual"):
                    if last_value < ax.get_ylim()[0]:
                        ax.set_ylim(ax.get_ylim()[0]*.1, ax.get_ylim()[1])
                    elif last_value > ax.get_ylim()[1]:
                        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*10)
                else:
                    if self._data.shape[0] == 2:
                        ax.set_ylim(self._data[col].min(), self._data[col].max())
                    else:
                        ylim = ax.get_ylim()
                        if last_value < ylim[0]:
                            ax.set_ylim(ylim[0] - .5*(ylim[1] - ylim[0]), ylim[1])
                        elif last_value > ylim[1]:
                            ax.set_ylim(ylim[0], ylim[1] + .5*(ylim[1] - ylim[0]))
        for i_col in range(len(self._plot_columns)):
            self._curves[i_col].set_data(self._data["iteration"], self._data[self._plot_columns[i_col]])
        return self._curves

class Basis:
    r"""! \brief Computes numerical parameters for nodal polynmial bases (such as `hexed::Gauss_legendre` and `hexed::Gauss_lobatto`)
    that are based on [Gaussian quadrature rules](https://en.wikipedia.org/wiki/Gaussian_quadrature).
    \details Used in `auto_generate.py`.
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
        bounds = np.zeros((self.row_size, 2))
        for i_side in range(2):
            for i_node in range(self.row_size):
                bounds[i_node, i_side] = self.interpolate(i_node, i_side)/self.weights[i_node]
        #bounds -= np.array(self.weights).astype(np.float64)@bounds
        def inner(vec0, vec1):
            return vec0@(np.array(self.weights).astype(np.float64)*vec1)
        self.space = np.zeros((self.row_size, 2))
        for i_side in range(2): bounds[:, i_side] /= np.sqrt(inner(bounds[:, i_side], bounds[:, i_side]))
        for i_side in range(2):
            self.space[:, i_side] = bounds[:, 0] + (1 - 2*i_side)*bounds[:, 1]
            self.space[:, i_side] /= np.sqrt(inner(self.space[:, i_side], self.space[:, i_side]))

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

    def get_ortho(self, degree, i_node, calc = False):
        r"""! equivalent to `self.legendre(degree)[i_node]` """
        return sp.Float(self.ortho[degree][i_node], self.calc_digits if calc else self.repr_digits)

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

    def filter(self, i_result, i_operand):
        dot = sp.Float(0, self.calc_digits)
        for i_inner in range(self.row_size):
            dot += self.get_ortho(i_inner, i_result, True)*0.5**i_inner*self.get_ortho(i_inner, i_operand, True)*self.weights[i_operand]
        return sp.Float(dot, self.repr_digits)

    def discretizations(self, n_elem):
        r"""! \brief Computes discrete operators for linear stability analsysis.
        \details Returns two matrices which are respectively the discrete analogues of
        \f$ -\partial u/\partial x \f$ and \f$ \partial^2 u/\partial x^2 \f$
        on a periodic 1D mesh with `n_elem` elements.
        """
        # sorry for the lack of comments... remind me to get back to this later
        nodes, weights = gauss_lobatto(self.row_size, self.calc_digits)
        nodes = [(node + 1)/2 for node in nodes]
        weights = [weight/2 for weight in weights]
        lobatto = Basis(nodes, weights, self.repr_digits, self.calc_digits)
        global_weights = np.zeros(n_elem*self.row_size)
        local_grad = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_avrg = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        neighb_jump = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
        discon = np.zeros((n_elem*self.row_size, n_elem*self.row_size))
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
                            discon[row, col] += (.5 - (i_side == j_side))*self.interpolate(i_col, 1-j_side)*lobatto.interpolate((1 - i_side)*(self.row_size - 1), self.nodes[i_row])
        advection = -local_grad + -neighb_avrg + neighb_jump
        diffusion = (local_grad + neighb_avrg)@(local_grad + neighb_avrg)
        assert np.linalg.norm(global_weights@advection) < 1e-12
        assert np.linalg.norm(global_weights@diffusion) < 1e-12
        return advection, diffusion

    def eigenvals(self, n_elem = 16):
        r"""! \brief Computes eigenvalues for use in time step calculation.
        \details Returns a list of two values which are the minimum real parts of the convection and diffusion operators, respectively.
        Eigenvalues are evaluated on a 1D, uniformly spaced, periodic mesh with `n_elem` elements for the linear advection/diffusion equations.
        These eigenvalues are nondimensionalized, so they apply to unit mesh spacing, wave speed, and diffusivity.
        \note standard floating-point precision (whatever that is for Python -- I think double)
        """
        return [np.linalg.eigvals(mat).real.min() for mat in self.discretizations(n_elem)]
