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

## \namespace hexed_python.std \brief wrapper for C++ standard library namespace `std`.

for path in include_paths:
    cppyy.add_include_path(path)
for lib in libraries:
    cppyy.load_library(lib)
cppyy.include("hexed/math.hpp")
cppyy.include("hexed/Solver_interface.hpp")
cppyy.include("hexed/Simplex_geom.hpp")
## \cond
cpp = cppyy.gbl.hexed
std = cppyy.gbl.std
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
    r"""! \brief Creates a `hexed::Solver_interface` object with a tree mesh.
    \details Mesh domain will be an axis-aligned box, optionally with pieces cut out of it by some user-specified geometry(s).
    Geometries need not be watertight.
    If any combination of the geometries creates a closed manifold with no gaps substantially larger than an element
    (in the fully refined mesh), then either the inside or the outside will be meshed, but not both
    (or possibly neither if there are nested closed geometries.
    For now, all sides of the box must have the same length.

    __Specifying resolution__

    Mesh resolution requirements may be specified as one of two ways:
    - As an `int`, which is interpreted as a lower bound on the refinement level
    - As a `float`, which is interpreted as an upper bound on the element side length

    \param n_dim (int) Number of dimenstions must be .
    \param row_size (int) Size of each row of quadrature points (total will be `row_size**n_dim` per element)
                    Must satisfy `2 <= row_size <= hexed::config::max_row_size`
    \param min_corner (float-array-like) minimum corner of the mesh bounding box. Must have size `n_dim`.
    \param max_corner (float-array-like) maximum corner of the mesh bounding box. Must have size `n_dim`.
    \param init_resolution (int or float) mesh will be refined to this resolution before geometry is inserted,
                           which usually also ends up being the farfield resolution.
                           If this is so coarse that inserting the geometry causes all elements to be deleted,
                           this will result in an exception.
                           The default value is usually acceptable for external aerodynamics applications.
    \param geometries (iterable) List of surface geometries to mesh, in one of the following formats:
                      - `n*2` numpy array representing the nodes of a polygonal curve (2D only)
    \param flood_fill_start (float-array-like or `None`) Seed point for flood fill algorithm.
                            That is, if any geometry in `geometries` divides the domain into disjoint regions, the region containing this point will be meshed.
                            Must have size `n_dim`.
                            If `None`, defaults to just inside `min_corner` (specifically `min_corner + 1e-6*(max_corner - min_corner)`)
    \param n_smooth (int) number of mesh smoothing iterations (meaning vertex movement, not refinement level smoothing)
                    to run after each refinement sweep (after the geometry is inserted).
    \param refine_sweeps Number of geometry-based refinement sweeps to run after geometry is inserted (_after_ `init_resolution` has already been achieved).
    \param surf_rep_tol Tolerance of surface representation (has to do with the discontinuity of surface normals between elements).
                        Smaller number means finer mesh.
                        To make this number completely irrelevant, specify something like 1000.
    \param surf_resolution Requirement on the resolution of elements at the geometry surface (usually leave as 0).
    \param max_resolution Elements will not be be refined past this resolution. The default value is fine enough that it effectively imposes no limit.
    """
    if not 1 <= int(n_dim) <= 3: raise User_error("invalid `n_dim`")
    if not 2 <= int(row_size) <= cpp.config.max_row_size: raise User_error(f"invalid `row_size` (max is {cpp.config.max_row_size})")
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
    #cpp_geoms = std.vector[std.unique_ptr[cpp.Surface_geom]]
    cpp_geoms = []
    for geom in geometries:
        if isinstance(geom, np.ndarray):
            if len(geom.shape) != 2: raise User_error("geometry array must be 2D")
            if (geom.shape[1] == 2):
                print("foo")
                cpp_geoms.append(cpp.Simplex_geom2(cpp.segments(to_matrix(geom))))
            else:
                raise User_error("geometry array must have 2 columns")
        else:
            raise User_error(f"could not interpret {type(geom)} as one of the supported geometry formats")
    for i_ref in range(ref_level(init_resolution)):
        solver.mesh().update()
    if cpp_geoms:
        if flood_fill_start is None:
            flood_fill_start = min_corner + 1e-6*(max_corner - min_corner)
        solver.mesh().set_surface(cpp.new_copy(cpp_geoms[0]), cpp.new_copy(cpp.Nonpenetration()), to_matrix(flood_fill_start))
    for i_smooth in range(n_smooth):
        solver.mesh().relax()
    solver.calc_jacobian()
    return solver

def naca(desig, n_points = 1000, closure = "warp"):
    r"""! \brief Constructs a NACA 4-digit airfoil geometry.
    \details Returns an n by 2 numpy array representing the coordinates of the airfoil at discrete points.
    This array can then be passed to `Solver.generate_mesh` as a geometry.
    Points are clustered near the leading edge but not the trailing (see implementation for details).
    This function is the recommended way to generate NACA airfoil geometry for Hexed simulations,
    as importing airfoils from coordinate files requires some \ref geom_fitting "special care".
    \param desig String representing the airfoil designation (e.g., "0012" for the NACA0012).
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
