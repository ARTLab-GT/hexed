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
    This is a family of airfoils with a cusped trailing edge which have analytic solutions for the incompressible flow around them.
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
    def infinite_generator(self):
        iteration = 0
        while not self.stop:
            yield iteration
            iteration += 1

    def get_new_lines(self):
        if not os.path.exists(self.directory + "output.txt"): return
        with open(self.directory + "output.txt", "rb") as output_file:
            output_file.seek(self.file_position, 0)
            self.lines += [line.decode("utf-8") for line in output_file.readlines()]
            self.file_position = output_file.tell()

    def __init__(self, directory = "hexed_out", interval = 0.2):
        if directory[-1] != "/": directory += "/"
        self.directory = directory
        self.interval = interval
        self.file_position = 0
        self.data = None
        self.lines = []
        while self.data is None:
            self.get_new_lines()
            while self.lines:
                line = self.lines.pop(0).replace(" ", "")
                if re.match("iteration,", line):
                    self.data = pd.DataFrame(columns = line.split(","))
                    break
            time.sleep(self.interval)
        self.fig = plt.figure()
        self.line, = plt.plot([], [])
        self.stop = False
        ani = FuncAnimation(self.fig, self.update, frames = self.infinite_generator, init_func = self.init,
                            blit = True, repeat = False, interval = int(self.interval*1e3), cache_frame_data = False)
        plt.show()

    def init(self):
        plt.xlim(0, 1)
        plt.ylim(0.1, 1.)
        plt.yscale("log")
        return self.line,

    def update(self, iteration):
        self.get_new_lines()
        while self.lines:
            line = self.lines.pop(0)
            if re.match(" *[0-9]+,", line):
                entries = line.split(",")
                self.data.loc[self.data.shape[0]] = [int(entries[0])] + [float(e) for e in entries[1:]]
                last_iter = self.data["iteration"][self.data.shape[0] - 1]
                last_value = self.data["normalized_residual"][self.data.shape[0] - 1]
                if last_iter > plt.xlim()[1]:
                    plt.xlim(0, plt.xlim()[1]*2)
                if last_value < plt.ylim()[0]:
                    plt.ylim(plt.ylim()[0]*.1, plt.ylim()[1])
                elif last_value > plt.ylim()[1]:
                    plt.ylim(plt.ylim()[0], plt.ylim()[1]*10)
        self.line.set_data(self.data["iteration"], self.data["normalized_residual"])
        return self.line,
