import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.optimize import fsolve
import os
import time
import re
import pandas as pd
import warnings

## \namespace hexedpy.utils
# \brief A Python module with random tools that are useful for Hexed simulations.

def naca(desig, n_points = 1000, closure = "warp"):
    r"""! \brief Constructs a NACA 4-digit airfoil geometry.
    \details Returns an n by 2 numpy array representing the coordinates of the airfoil at discrete points.
    This array can then be passed to `Solver.generate_mesh` as a geometry.
    Points are clustered near the leading edge but not the trailing (see implementation for details).
    This function is the recommended way to generate NACA airfoil geometry for Hexed simulations,
    as importing airfoils from coordinate files requires some \ref geom_fitting "special care".
    \param desig String representing the airfoil designation (e.g., `"0012"` for the NACA0012).
        In general, we cannot accept this parameter as an `int` because of possible leading zeros.
    \param n_points Number of points on the airfoil surface. Don't be stingy,
        since DG is finnicky with discrete geometry representations---1000 is actually on the lower end of what I normally use.
    \param closure If and how to close the trailing edge. There are 3 options:
        - `"warp"`: Close the trailing edge by adding a 4th-degree polyomial of \f$ x_0 \f$ to \f$ x_1 \f$.
        - `"segment"`: Close the trailing edge by adding a line segment connecting the last point to the first point,
            causing the array to be `(n_points + 1)*2` instead of `n_points*2`.
        - `"extend"`: Extends the \f$ x_0 \f$ domain until the point where the upper and lower surface meet
            and then rescales the whole airfoil to keep the chord equal to 1.
            This stays faithful to the original polynomial profile, but slightly decreases the thickness.
            This is the method used in the classic [NASA validation case](https://turbmodels.larc.nasa.gov/naca0012_val.html).
            Extension and rescaling are applied __after__ the camber addition.
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
    def sym_profile(p):
        ap = np.abs(p)
        return 5*thickness*p*(.2969 - .1260*ap - .3516*ap**3 + .2843*ap**5 - .1015*ap**7)
    max_param = 1.
    if (closure == "extend"):
        max_param = fsolve(sym_profile, 1., xtol=1e-12)[0]
    param = np.linspace(-max_param, max_param, n_points)
    coords[:, 0] = param**2
    coords[:, 1] = sym_profile(param)
    if camber_loc > 0:
        camber_max *= 1e-2
        camber_loc *= 1e-1
        s = coords[:, 0] <  camber_loc
        coords[s, 1] += camber_max/camber_loc**2*(2*camber_loc*param[s]**2 - param[s]**4)
        s = coords[:, 0] >= camber_loc
        coords[s, 1] += camber_max/(1 - camber_loc)**2*(1 - 2*camber_loc + 2*camber_loc*param[s]**2 - param[s]**4)
    if closure == "warp":
        coords[:, 1] -= param*np.abs(param)**7*coords[-1, 1]
    elif closure == "segment":
        coords = np.concatenate([coords, coords[[0], :]])
    elif closure == "extend":
        coords /= max_param**2
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

def sort_curve(data, coords, start_at=0, closed=True):
    r"""! \brief Returns data sorted into a continuous curve.
    \param data A [pandas.DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html)
        containing the data (coordinates and flow variables) to be sorted
    \param coords Labels of columns to treat as coordinates
    \param start_at Index of point to start at
    \param closed If true, will append an additional copy of the start point to the end
        to ensure that the final curve is closed.
    """
    old_inds = list(range(len(data)))
    new_inds = [old_inds.pop(start_at)]
    while len(old_inds):
        def get_dist(ind):
            return np.linalg.norm(data[coords].iloc[ind] - data[coords].iloc[new_inds[-1]])
        nearest = np.argmin([get_dist(ind) for ind in old_inds])
        new_inds.append(old_inds.pop(nearest))
    if closed:
        new_inds.append(new_inds[0])
    return data.reindex(index=new_inds)

def read_history(output_directory):
    r"""! \brief Reads the convergence history from `output.txt` and returns it as a
    [pandas.DataFrame](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.html).
    \param output_directory The \ref working_dir of the simulation you want to read the history of
        (e.g. `hexed_out` if you used the default `working_dir`).
        Can be an absolute or relative path.
    """
    fname = f"{output_directory}/output.txt"
    with open(fname, "r") as output_file:
        text = output_file.read()
    lines = text.split("\n")
    while lines and "Meshing complete" not in lines[0]:
        lines.pop(0)
    while lines and "iteration" not in lines[0]:
        lines.pop(0)
    assert lines, f"Could not find any convergence history data in `{fname}`."
    col_names = [name.strip() for name in lines.pop(0).split(",")]
    col_data = []
    for line in lines[1:]:
        if "simulation complete" in line:
            break
        if line.startswith(" "):
            row = []
            for entry in line.split(","):
                if all([c.isnumeric() or c.isspace() for c in entry]):
                    row.append(int(entry))
                else:
                    row.append(float(entry))
        col_data.append(row)
    return pd.DataFrame(data = col_data, columns = col_names)

class History_plot:
    r"""! \brief creates a real-time, interactive plot of the convergence history
    \details Convergence history is obtained from the `output.txt` file which contains the console output of \ref hexecute.
    Every column in the output whose name is not in `History_plot.column_blacklist` will be plotted in its own subplot.
    Columns with names ending in `residual` or `error` will be plotted on a log scale.
    A `History_plot` instance should be created in a separate process from the solver
    but may be created before, after, or during the simulation---the plot will not appear until an `output.txt` file exists.
    A `History_plot` can be created directly from the \ref hil "HIL" solver script with $\ref plot_history.
    """
    ## \brief names of output columns __not__ to plot
    ## \details You may modify this variable for the class or for instances
    column_blacklist = ["flow_time", "time_step", "time_stage"]

    def _infinite_generator(self):
        while not self._stop:
            yield None

    def _read_new_lines(self):
        if not os.path.exists(self._directory + "output.txt"): return
        with open(self._directory + "output.txt", "rb") as output_file:
            output_file.seek(self._file_position, 0)
            self._lines += [line.decode("utf-8") for line in output_file.readlines()]
            self._file_position = output_file.tell()

    def __init__(self, directory = "hexed_out", interval = 2.):
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
        self._monitor_window = None
        while self._data is None:
            self._read_new_lines()
            while self._lines:
                line = self._lines.pop(0).replace(" ", "").replace("\n", "")
                if line.startswith("Monitorwindow:"):
                    self._monitor_window = float(line.split(":")[1])
                if re.match("iteration,", line):
                    self._data = pd.DataFrame(columns = line.split(","))
                    break
            time.sleep(self._interval)
        assert self._monitor_window is not None, "Monitor window is not indicated in output text."
        def plot_column(col):
            if col in self.column_blacklist:
                return False
            if "iteration" in col:
                return False
            return True
        self._plot_columns = [col for col in self._data.columns if plot_column(col)]
        self._stop = False
        self._fig, self._axs = plt.subplots(1, len(self._plot_columns))
        plt.tight_layout()
        self._fig.set_size_inches(18, 5)
        ani = FuncAnimation(self._fig, self._update, frames = self._infinite_generator, init_func = self._init,
                            blit = False, repeat = False, interval = int(self._interval*1e3), cache_frame_data = False)
        plt.show()

    def _init(self):
        self._curves = []
        self._stats = {}
        for i_col in range(len(self._plot_columns)):
            col = self._plot_columns[i_col]
            label = col.replace("_", " ")
            if label == "pseudotime iteration":
                label = label + "s"
            ax = self._axs[i_col]
            self._curves.append(ax.plot([], [])[0])
            ax.set_xlim(0., 1.)
            ax.grid(True)
            ax.set_xlabel("iteration")
            ax.set_ylabel(label)
            if label.endswith("residual") or label.endswith("error"):
                self._axs[i_col].set_ylim(0.1, 1.)
                self._axs[i_col].set_yscale("log")
            self._stats[col] = [
                ax.plot([], [], color="black")[0],
                ax.plot([], [], color="grey", linestyle="dashed")[0],
                ax.plot([], [], color="grey", linestyle="dashed")[0],
            ]
        return self._curves

    def _update(self, _):
        self._read_new_lines()
        while self._lines:
            line = self._lines.pop(0)
            if line.startswith("simulation complete"):
                self._stop = True
            try:
                status_data = pd.read_csv(self._directory + "status_data.txt", delimiter=":", names=["parameter", "value"], index_col=0)
                has_status = True
            except:
                has_status = False
            if re.match(" *[0-9]+,", line):
                entries = line.split(",")
                add_line = self._data.shape[0]
                if add_line > 0 and int(entries[0]) == self._data["iteration"][add_line - 1]:
                    add_line -= 1
                self._data.loc[add_line] = [int(entries[0])] + [float(e) for e in entries[1:]]
                last_iter = self._data["iteration"][self._data.shape[0] - 1]
                if last_iter > self._axs[0].get_xlim()[1]:
                    for ax in self._axs:
                        ax.set_xlim(0, self._data["iteration"].max()*2)
                for i_col in range(len(self._plot_columns)):
                    ax = self._axs[i_col]
                    col = self._plot_columns[i_col]
                    last_value = self._data[col][self._data.shape[0] - 1]
                    if col.endswith("residual") or col.endswith("error"):
                        if last_value < ax.get_ylim()[0]:
                            ax.set_ylim(self._data[col].min()*.1, ax.get_ylim()[1])
                        elif last_value > ax.get_ylim()[1]:
                            ax.set_ylim(ax.get_ylim()[0], self._data[col].max()*10)
                    else:
                        if self._data.shape[0] == 2:
                            ax.set_ylim(self._data[col].min(), self._data[col].max())
                        else:
                            ylim = ax.get_ylim()
                            if last_value < ylim[0]:
                                ax.set_ylim(ylim[1] + 1.5*(self._data[col].min() - ylim[1]), ylim[1])
                            elif last_value > ylim[1]:
                                ax.set_ylim(ylim[0], ylim[0] + 1.5*(self._data[col].max() - ylim[0]))
                    if has_status:
                        if col + "_smoothed" in status_data.index:
                            smoothed = status_data.at[col + "_smoothed", "value"]
                            trend = status_data.at[col + "_trend", "value"]
                            noise = status_data.at[col + "_noise", "value"]
                            noise_trend = status_data.at[col + "_noise_trend", "value"]
                            iteration = self._data.at[add_line, "iteration"]
                            x = [(1 - self._monitor_window)*iteration, iteration]
                            y = np.array([smoothed - trend*self._monitor_window*iteration, smoothed]);
                            self._stats[col][0].set_data(x, y)
                            spread = np.array([noise - noise_trend*self._monitor_window*iteration, noise])
                            self._stats[col][1].set_data(x, y - spread)
                            self._stats[col][2].set_data(x, y + spread)
        for i_col in range(len(self._plot_columns)):
            self._curves[i_col].set_data(self._data["iteration"], self._data[self._plot_columns[i_col]])
        return self._curves
