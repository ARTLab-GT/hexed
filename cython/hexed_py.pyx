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

    def generate_mesh(self, min_corner, max_corner, init_resolution = 3):
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

        \param min_corner (float-array-like) minimum corner of the mesh bounding box. Must have size `n_dim`
        \param max_corner (float-array-like) maximum corner of the mesh bounding box. Must have size `n_dim`
        \param init_resolution (int or float) mesh will be refined to this resolution before geometry is inserted,
                               which usually also ends up being the farfield resolution.
                               If this is so coarse that inserting the geometry causes all elements to be deleted,
                               this will result in an exception.
                               The default value is usually acceptable for external aerodynamics applications.
        """
        if self._is_init:
            raise User_error("this `Solver` already has a mesh")
        try:
            min_corner = np.array(min_corner).astype(np.float64).flatten()
            max_corner = np.array(max_corner).astype(np.float64).flatten()
        except:
            raise User_error("could not convert `min_corner` and `max_corner` to float arrays")
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
        if isinstance(init_resolution, int):
            ref_level = init_resolution
        elif isinstance(init_resolution, float):
            ref_level = math.ceil(np.log(root_sz/init_resolution)/np.log(2))
        else:
            raise User_error("could not interpret `init_resolution` as `int` or `float`")
        ref_level = max(ref_level, 0)
        for i_ref in range(ref_level):
            self._solver[0].mesh().update()

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
        self._solver[0].calc_jacobian()
        self._solver[0].visualize_field_tecplot(bytes(f"{self._output_dir}/{file_name}", "ascii"), n_sample, False, False, True)
