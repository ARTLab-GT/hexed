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

def to_arr(arr_like, shape = None, size = None, exception = User_error("Could not convert to numpy float array")):
    r"""! converts `arr_like` to a numpy float array, potentially with a specific size and/or
    shape and throws a specific exception on failure """
    try:
        arr = np.array(arr_like).astype(np.float64)
        if shape is not None: assert arr.shape == shape
        if size is not None: assert arr.size == size
        return arr
    except Exception as e:
        raise exception from e

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

## \cond
class _Criterion:
    instances = 0
## \endcond

def ref_criterion(requirement):
    r"""! \brief Obtains a C++ function that can be used to determine which elements should be refined.
    \details `requirement` is something that can somehow be interpreted as a formula to look at an element and get a boolean.
    The following data types are understood:
    - `bool`: Gives this value for every element.
    - `int`: Refinement level bound -- gives `True` for any element with refinement level __less than__ this.
    - `float`: Size bound -- gives `True` for any element with side length __greater than__ this
    - `str`: A C++ expression which will be evaluated for each element in an environment containing the following variables:
      - __bool is_extruded:__ Whether the element is extruded, which basically means whether it is on the geometry surface.
      - __int ref_level:__ Current refinement level of the element.
      - __double nom_sz:__ Current nominal size of the element (the size the element _would_ be if it weren't deformed).
      - __Eigen::Vector<double, 3, 1> center:__ A vector object representing position of the center of the element
        (technically, the average position of its vertices).
        This vector has size 3 (trailing dimensions set to 0), supports basic arithmetic operators, and has `.dot(other)` and `.norm()` methods.
      - __double res_bad:__ The `hexed::Element::resolution_badness` parameter.
      You also have access to all the usual C++ featurs including standard library math functions (e.g. `std::exp`, `std::log`)
      and the logical operators `&&` and `||`.
      Under the hood, your expression will be the body of an anonymous function
      compiled Just In Time with cling and passed to `hexed::criteria::criterion`.
    """
    if requirement is False:
        return cpp.criteria.never
    elif requirement is True:
        return cpp.criteria.always
    elif isinstance(requirement, int):
        return ref_criterion(f"ref_level < {requirement}")
    elif isinstance(requirement, float):
        return ref_criterion(f"nom_sz > {requirement}")
    elif isinstance(requirement, str):
        name = f"criterion{_Criterion.instances}"
        _Criterion.instances += 1
        cppyy.cppdef(f"""namespace hexed::criteria {{
        auto {name} = hexed::criteria::criterion(
            [](bool is_extruded, int ref_level, double nom_sz, Eigen::Vector3d center, double res_bad) {{
                return bool({requirement});
        }});
        }}""")
        return getattr(cpp.criteria, name)
    else:
        raise User_error(f"invalid refinement criterion type {type(requirement)}")

def create_solver(
        n_dim, row_size,
        min_corner, max_corner,
        init_cond,
        geometries = [], flood_fill_start = None, n_smooth = 10,
        init_resolution = 3, init_max_iters = 1000,
        final_resolution = False, final_max_iters = 1000,
    ):
    r"""! \brief Creates a `hexed::Solver_interface` object with a premade tree mesh.
    \details This is the recommended way for end users to construct a Solver object.
    You can still interact with the `mesh()` method to further adjust the mesh, but you shouldn't need to.
    Mesh domain will be an axis-aligned box, optionally with pieces cut out of it by some user-specified geometry(s).
    Geometries need not be watertight.
    If any combination of the geometries creates a closed manifold with no gaps substantially larger than an element
    (in the fully refined mesh), then either the inside or the outside will be meshed, but not both
    (or possibly neither if there are nested closed geometries.
    For now, all sides of the box must have the same length.

    \param n_dim (int) Number of dimenstions must be .
    \param row_size (int) Size of each row of quadrature points (total will be `row_size**n_dim` per element)
                    Must satisfy `2 <= row_size <= hexed::config::max_row_size`
    \param min_corner (array-like) minimum corner of the mesh bounding box. Must have size `n_dim`.
    \param max_corner (array-like) maximum corner of the mesh bounding box. Must have size `n_dim`.
    \param init_cond (hexed::Spacetime_func or array-like) Specifies the state variables at time 0
                     in the \ref state_vector "momentum-density-energy order". If an array is passed it will be converted to a `hexed::Constant_func`.
    \param geometries (iterable) List of surface geometries to mesh, in one of the following formats:
                      - `n*2` numpy array representing the nodes of a polygonal curve (2D only)
    \param flood_fill_start (array-like or `None`) Seed point for flood fill algorithm.
                            That is, if any geometry in `geometries` divides the domain into disjoint regions, the region containing this point will be meshed.
                            Must have size `n_dim`.
                            If `None`, defaults to just inside `min_corner` (specifically `min_corner + 1e-6*(max_corner - min_corner)`)
    \param n_smooth (int) number of mesh smoothing iterations (meaning vertex movement, not refinement level smoothing)
                    to run after each refinement sweep (after the geometry is inserted).
    \param init_resolution Argument to `ref_criterion()` to control the resolution before the geometry is inserted,
                           which usually also ends up being the farfield resolution.
                           The mesh will be iteratively refined until no more elements want to refine or `init_max_iters` is reached.
                           If this is so coarse that inserting the geometry causes all elements to be deleted,
                           this will result in an exception.
                           The default value is usually acceptable for external aerodynamics applications.
    \param init_max_iters (int) Maximum number of refinement iterations to execute before geometry is inserted.
    \param final_resolution Argument to `ref_criterion()` to control the resolution after the geometry is inserted,
                            The mesh will be iteratively refined until no more elements want to refine or `final_max_iters` is reached.
                            Before each refinement iteration, the `hexed::Element::resolution_badness` will be set to a value
                            which reflects the quality of the surface representation (large number means surface is poorly resolved).
    \param final_max_iters (int) Maximum number of refinement iterations to execute after geometry is inserted.
    """
    if not 1 <= int(n_dim) <= 3: raise User_error("invalid `n_dim`")
    if not 2 <= int(row_size) <= cpp.config.max_row_size: raise User_error(f"invalid `row_size` (max is {cpp.config.max_row_size})")
    min_corner = to_arr(min_corner, size = n_dim).flatten()
    max_corner = to_arr(max_corner, size = n_dim).flatten()
    if not np.all(np.array(min_corner) < np.array(max_corner)): raise User_error(f"`min_corner` must be strictly less than `max_corner`")
    root_sz = np.max(max_corner - min_corner)
    solver = cpp.make_solver(n_dim, row_size, root_sz)
    solver.mesh().add_tree([cpp.new_copy(cpp.Nonpenetration()) for i in range(2*n_dim)], to_matrix(min_corner))
    crit = ref_criterion(init_resolution)
    for i_refine in range(init_max_iters):
        if not solver.mesh().update(crit):
            break
    cpp_geoms = []
    for geom in geometries:
        if isinstance(geom, np.ndarray):
            if len(geom.shape) != 2: raise User_error("geometry array must be 2D")
            if (geom.shape[1] == 2):
                cpp_geoms.append(cpp.Simplex_geom2(cpp.segments(to_matrix(geom))))
            else:
                raise User_error("geometry array must have 2 columns")
        else:
            raise User_error(f"could not interpret {type(geom)} as one of the supported geometry formats")
    if cpp_geoms:
        if flood_fill_start is None:
            flood_fill_start = min_corner + 1e-6*(max_corner - min_corner)
        solver.mesh().set_surface(cpp.new_copy(cpp_geoms[0]), cpp.new_copy(cpp.Nonpenetration()), to_matrix(flood_fill_start))
    for i_smooth in range(n_smooth):
        solver.mesh().relax()
    if solver.mesh().n_elements() == 0: raise User_error("geometry addition deleted all elements")
    crit = ref_criterion(final_resolution)
    for i_refine in range(final_max_iters):
        if not solver.mesh().update(crit):
            break
        for i_smooth in range(n_smooth):
            solver.mesh().relax()
    solver.calc_jacobian()
    if not isinstance(init_cond, cpp.Spacetime_func):
        init_cond = to_arr(init_cond, size = n_dim + 2, exception = User_error(f"could not interpret `init_cond` as `Spacetime_func` or array of size {n_dim + 2}"))
        init_cond = cpp.Constant_func(init_cond.flatten())
    solver.initialize(init_cond);
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

## \cond
# altitude vs temperature
_standard_atm_temp = np.array([
[    0., 288.15],
[11000., 216.65],
[20000., 216.65],
[32000., 228.65],
[47000., 270.65],
[51000., 270.65],
[71000., 214.65],
[80000., 196.65],
])
## \endcond

def standard_atm(alt_geom, temp_offset = 0.):
    r"""! \brief computes the [ICAO Standard atmosphere](http://www.aviationchief.com/uploads/9/2/0/9/92098238/icao_doc_7488_-_manual_of_icao_standard_atmosphere_-_3rd_edition_-_1994.pdf)
    \details Valid from 0 to 80km.
    \param alt_geom Geometric altitude (as opposed to geopotential altitude)
    \param temp_offset Temperature will be incremented by `temp_offset` relative to the standard atmosphere without changing the pressure.
    \returns pair (density, pressure)
    \todo is this temperature offset behavior correct?
    \see \ref units
    """
    # convert geometric altitude to geopotential altitude
    alt_geopot = cpp.earth_radius*alt_geom/(cpp.earth_radius + alt_geom)
    pres = cpp.atmosphere
    n_rows = (_standard_atm_temp[:, 0] <= alt_geopot + 1e-12).sum() # small epsilon added to prevent errors at ~0 altitude
    if n_rows > _standard_atm_temp.shape[0]: raise User_error("altitude is above model limit")
    elif n_rows == 0: raise User_error("altitude is below model limit")
    for i_row in range(n_rows):
        base_temp = _standard_atm_temp[i_row, 1]
        lapse = (_standard_atm_temp[i_row + 1, 1] - base_temp)/(_standard_atm_temp[i_row + 1, 0] - _standard_atm_temp[i_row, 0])
        if i_row == n_rows - 1:
            h = alt_geopot
        else:
            h = _standard_atm_temp[i_row + 1, 0]
        h_diff = h - _standard_atm_temp[i_row, 0]
        temp = base_temp + lapse*h_diff
        if abs(lapse) < 1e-10:
            pres *= np.exp(-cpp.std_grav*h_diff/cpp.specific_gas_air/base_temp)
        else:
            pres *= (temp/base_temp)**(-cpp.std_grav/lapse/cpp.specific_gas_air)
    temp += temp_offset
    dens = pres/cpp.specific_gas_air/temp
    return dens, pres

def rot_mat(angle):
    r"""! \brief compute a 2D rotation matrix """
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

def flow_state(density = None, pressure = None, temperature = None, altitude = None,
               mach = None, speed = None,
               direction = None, velocity = None, attack = None, sideslip = None):
    r"""! \brief Computes the \ref state_vector based on parameters of practical engineering relevance.
    \details Specify some combination of the keyword arguments that fully defines the state.
    The dimensionality is inferred from how you specify the direction:
    - If you specify (angle of) `attack` but not `sideslip`, you get a 2D state.
      If you specify both, you get 3D.
    - If you specify a `direction` or `velocity`, the number of components determines the dimensionality.
    If the state is underdetermined you will get an exception.
    If overdetermined, it will choose which constraints to satisfy.
    Assumes the flow is CPG air with heat ratio \f$ \gamma = 1.4 \f$.

    \attention Because all quantities are in SI base units, the angles are radian!
    To specify angles in degrees, you can do it like:
    ```
    attack = 10*hexed.cpp.degree
    ```
    \see \ref units
    \see `constants.hpp`
    """
    heat_ratio = 1.4
    if altitude is not None:
        density, pressure = standard_atm(altitude)
    try:
        # find the density and pressure
        # note that the arithmetic operations will throw exceptions if any of their arguments are None
        if density is None:
            density = pressure/cpp.specific_gas_air/temperature
        elif pressure is None:
            pressure = density*cpp.specific_gas_air*temperature
        int_ener = pressure/(heat_ratio - 1.)
        # find flow direction, if applicable
        if attack is not None:
            if sideslip is None:
                direction = rot_mat(attack)@np.array([1, 0])
            else:
                direction = np.array([1, 0, 0])
                direction[1:] = rot_mat(-attack) @ direction[1:]
                direction[:-1] = -rot_mat(sideslip) @ direction[:-1]
        # find the velocity
        if velocity is None:
            if mach is not None:
                speed = mach*(heat_ratio*pressure/density)**.5
            velocity = speed*direction
    except Exception as e:
        raise User_error("underdetermined flow state specification") from e
    return np.concatenate([density*velocity, [density, pressure/(heat_ratio - 1.) + .5*density*velocity@velocity]])
