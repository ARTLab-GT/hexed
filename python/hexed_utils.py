import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import time
import re
import pandas as pd

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
                line = self._lines.pop(0).replace(" ", "")
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
            ax = self._axs[i_col]
            self._curves.append(ax.plot([], [])[0])
            ax.set_xlim(0, 1)
            ax.grid(True)
            ax.set_xlabel("iteration")
            ax.set_ylabel(self._plot_columns[i_col])
            if self._plot_columns[i_col].endswith("residual"):
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
