import cppyy
import numpy as np
import math
from hexed_py_config import *

## \namespace hexed_python
# \brief Internal namespace for the Python API.
# \details In the Doxygen documentation,
# everything in the Python API will appear in the namespace `hexed_python` because it is defined in the file `python/hexed_python.py`
# However, when you actually use it, it lives in the module `hexed` (import it with `import hexed`).
# This distinction prevents naming ambiguities in the documentation,
# since both the C++ namespace and the Python module are called `hexed`.

## \namespace hexed_python.cpp
# \brief Namespace for accessing the C++ API in Python.
# \details This submodule provides Python bindings to the C++ namespace `hexed`
# via [cppyy](https://cppyy.readthedocs.io/en/latest/starting.html).
# For the most part, everything should work pretty much exactly as it would in C++.
# For example, to call `hexed::math::pow``(3, 2)` in Python you can do `hexed.cpp.math.pow(3, 2)`.
# However, please read the docs on cppyy if you want to do anything fancy.
# This submodule doesn't provide access to all the C++ definitions,
# but it should provide everything you need for high-level solver execution.

for path in include_paths:
    cppyy.add_include_path(path)
for lib in libraries:
    cppyy.load_library(lib)
cppyy.include("hexed/math.hpp")
cppyy.include("hexed/Solver_interface.hpp")
## \cond
cpp = cppyy.gbl.hexed
## \endcond

class User_error(Exception):
    r"""! \brief An exception indicating that some user input/action was demonstrably invalid.
    \see \ref user_errors
    """
    pass

def to_arr(arr_like, shape = None, size = None, exception = User_error):
    r"""! converts `arr_like` to a numpy float array, potentially with a specific size and/or
    shape and throws a specific exception on failure """
    try:
        arr = np.array(arr_like).astype(np.float64)
        if shape is not None: assert arr.shape == shape
        if size is not None: assert arr.size == size
        return arr
    except Exception as e:
        raise exception("Could not convert to numpy float array") from e

def to_matrix(arr_like):
    r"""! Converts `arr_like` to an `Eigen::MatrixXd`. Array is transposed to maintain storage order. """
    arr = np.array(arr_like).astype(np.float64)
    if len(arr.shape) > 2:
        raise Exception("array has too many dimensions to be converted to a matrix")
    while len(arr.shape) < 2:
        arr = np.array([arr])
    mat = cppyy.gbl.Eigen.MatrixXd(arr.shape[1], arr.shape[0])
    view = mat.data()
    view.reshape((arr.size,))
    view[:] = arr.flatten()
    return mat

def create_solver(
        n_dim, row_size,
        min_corner, max_corner, init_resolution = 3,
        geometries = [], flood_fill_start = None, n_smooth = 10,
        refine_sweeps = 0, surf_rep_tol = 1000., surf_resolution = 0, max_resolution = 1000,
    ):
    min_corner = to_arr(min_corner, size = n_dim).flatten()
    max_corner = to_arr(max_corner, size = n_dim).flatten()
    if not np.all(np.array(min_corner) < np.array(max_corner)): raise User_error(f"`min_corner` must be strictly less than `max_corner`")
    root_sz = np.max(max_corner - min_corner)
    solver = cpp.make_solver(n_dim, row_size, root_sz)
    solver.mesh().add_tree([cpp.new_copy(cpp.Nonpenetration()) for i in range(2*n_dim)], to_matrix(min_corner))
    def ref_level(resolution):
        if isinstance(resolution, int):
            rl = resolution
        elif isinstance(resolution, float):
            rl = math.ceil(np.log(root_sz/resolution)/np.log(2))
        else:
            raise User_error(f"could not interpret resolution {resolution} as `int` or `float`")
        rl = max(rl, 0)
        return rl
    for i_ref in range(ref_level(init_resolution)):
        solver.mesh().update()
    return solver
