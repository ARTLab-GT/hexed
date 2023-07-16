import numpy as np
import math
import os
cimport hexed_cpp as cpp

## \namespace hexed_py
# \brief Internal namespace for the Python API.
# \details In the Doxygen documentation,
# everything in the Python API will appear in the namespace `hexed_py` because it is defined in the file `cython/hexed_py.pyx`.
# However, when you actually use it, it lives in the module `hexed` (import it with `import hexed`).
# This distinction prevents naming ambiguities in the documentation,
# since both the C++ namespace and the Python module are called `hexed`.

def to_np_arr(arr_like):
    r"""! converts `arr_like` to a numpy float array and throw a `User_error` on failure """
    try:
        return np.array(arr_like).astype(np.float64)
    except Exception as e:
        raise User_error("Could not convert to numpy float array") from e

def matrix_shape(arr):
    r"""! \brief Finds the shape a matrix must have to match the storage order of a <= 2D array """
    shape = list(arr.shape[::-1]) # size that matrix will have. reverse due to col vs row major discrepancy
    assert(len(shape) <= 2, "cannot convert a >2D array shape to a 2D matrix shape")
    # matrix must be 2D, so if array dimensionality is <2, set trailing dimensions to 1
    while len(shape) < 2:
        shape.append(1)
    return shape

cdef cpp.MatrixXd matrix(arr):
    arr = np.array(arr).astype(np.float64)
    shape = matrix_shape(arr)
    cdef cpp.MatrixXd mat = cpp.MatrixXd(shape[0], shape[1]);
    arr = arr.flatten()
    cdef double [:] arr_view = arr
    cdef double [:] mat_view = <double[:mat.size()]>mat.data()
    mat_view[:] = arr_view
    return mat

cdef class Iteration_status:
    r"""! \brief Python interface for `hexed::Iteration_status` """
    cdef cpp.Iteration_status status
    def header(self):
        return self.status.header().decode()
    def report(self):
        return self.status.report().decode()

class User_error(Exception):
    r"""! \brief An exception indicating that some user input/action was demonstrably invalid.
    \see \ref user_errors
    """
    pass

cdef class Solver:
    r"""! \brief Python interface for `hexed::Solver` """
    cdef cpp.Solver* _solver
    cdef bint _is_init # whether `_solver` has been initialized to point to some data
    cdef int _n_dim
    cdef int _row_size
    cdef int _lts
    cdef str _output_dir

    def __cinit__(self, int n_dim, int row_size, bint local_time_stepping = False, output_dir = "./hexed_out"):
        """!
        \param n_dim (int) Number of dimenstions. Must satisfy `1 <= n_dim <= 3`
        \param row_size (int) Size of each row of quadrature points (total will be `row_size**n_dim` per element)
                        Must satisfy `2 <= row_size <= hexed::config::max_row_size`
        \param local_time_stepping (bool) Whether solver should use local or global time stepping
        \param output_dir where to write output files. If this directory does not exist, it will be created
        """
        _is_init = False
        if not 1 <= int(n_dim) <= 3: raise User_error("invalid `n_dim`")
        if not 2 <= int(row_size) <= cpp.max_row_size: raise User_error(f"invalid `row_size` (max is {cpp.max_row_size})")
        self._n_dim = int(n_dim)
        self._row_size = int(row_size)
        self._lts = bool(local_time_stepping)
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        self._output_dir = output_dir
    def __dealloc__(self):
        if self._is_init:
            del self._solver

    def generate_mesh(self,
                      min_corner, max_corner, init_resolution = 3,
                      geometries = [], flood_fill_start = None, n_smooth = 10,
                      refine_sweeps = 0, surf_rep_tol = 1000., surf_resolution = 0, max_resolution = 1000,
                     ):
        r"""! \brief Generates a tree mesh for the simulation.
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
        if self._is_init:
            raise User_error("this `Solver` already has a mesh")
        min_corner = to_np_arr(min_corner)
        max_corner = to_np_arr(max_corner)
        if not min_corner.size == self._n_dim: raise User_error(f"`min_corner` must have size {self._n_dim}")
        if not max_corner.size == self._n_dim: raise User_error(f"`max_corner` must have size {self._n_dim}")
        if not np.all(np.array(min_corner) < np.array(max_corner)): raise User_error(f"`min_corner` must be strictly less than `max_corner`")
        root_sz = np.max(max_corner - min_corner)
        self._solver = new cpp.Solver(self._n_dim, self._row_size, root_sz, self._lts, cpp.inviscid, cpp.inviscid)
        self._is_init = True
        cdef cpp.vector[cpp.Flow_bc*] bcs
        for i_bc in range(2*self._n_dim):
            bcs.push_back(new cpp.Nonpenetration())
        self._solver[0].mesh().add_tree(bcs, matrix(min_corner))
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
            self._solver[0].mesh().update()
        cdef cpp.vector[cpp.unique_ptr[cpp.Surface_geom]] cpp_geoms
        for geom in geometries:
            if isinstance(geom, np.ndarray):
                if len(geom.shape) != 2: raise User_error("geometry array must be 2D")
                if (geom.shape[1] == 2):
                    cpp_geoms.push_back(cpp.unique_ptr[cpp.Surface_geom](new cpp.Simplex_geom2(cpp.segments(matrix(geom)))))
                else:
                    raise User_error("geometry array must have 2 columns")
            else:
                raise User_error(f"could not interpret following geometry specification as one of the supported formats:\n{geom}")
        if cpp_geoms.size() > 0:
            if flood_fill_start is None:
                flood_fill_start = min_corner + 1e-6*(max_corner - min_corner)
            self._solver[0].mesh().set_surface(cpp_geoms[0].release(), new cpp.Nonpenetration(), matrix(to_np_arr(flood_fill_start)))
        for i_smooth in range(n_smooth):
            self._solver[0].mesh().relax()
        for i_ref in range(refine_sweeps):
            self._solver[0].calc_jacobian()
            self._solver[0].set_res_bad_surface_rep(self._solver[0].mesh().surface_bc_sn())
            self._solver[0].mesh().update(cpp.Mesh.General_ref_criterion(surf_rep_tol, ref_level(surf_resolution), ref_level(max_resolution)))
        for i_smooth in range(n_smooth):
            self._solver[0].mesh().relax()
        self._solver[0].calc_jacobian()

    def iteration_status(self):
        r"""! \brief fetches the `Iteration_status` describing the state of the simulation """
        status = Iteration_status()
        status.status = self._solver[0].iteration_status()
        return status

    def visualize_field(self, file_name, int n_sample = 20):
        r"""! \brief Visualizes the field data.
        \details Data is written in Tecplot subzone data format (`.szplt`).
        \param file_name name of file to write the data to, without file type extension.
                         will be prefixed with `self.output_dir`
        \param n_sample each element will be a `n_sample`[*`n_sample`[*`n_sample]]
                        array of uniformly spaced sample points
        """
        self._solver[0].visualize_field_tecplot(bytes(f"{self._output_dir}/{file_name}", "ascii"), n_sample, False, False, True)

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
