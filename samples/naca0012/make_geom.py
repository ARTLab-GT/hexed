"""
Subsonic NACA 0012 Example Problem: Geometry generation script.
Creates a file `auto.csv` with airfoil coordinates.
This script is executed by the input file `run.hil`.
"""

import numpy as np
from hexedpy.utils import naca

data = naca("0012", 2**15 + 1)
np.savetxt("auto.csv", data, delimiter = ", ")
