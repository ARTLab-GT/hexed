import numpy as np
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

    def __cinit__(self, int n_dim, int row_size, bint local_time_stepping = False):
        """!
        \param n_dim (`int`) Number of dimenstions. Must satisfy `1 <= n_dim <= 3`
        \param row_size (`int`) Size of each row of quadrature points (total will be `row_size**n_dim` per element)
                        Must satisfy `2 <= row_size <= hexed::config::max_row_size`
        \param local_time_stepping (`bool`) Whether solver should use local or global time stepping
        """
        _is_init = False
        if not 1 <= int(n_dim) <= 3: raise User_error("invalid `n_dim`")
        if not 2 <= int(row_size) <= cpp.max_row_size: raise User_error(f"invalid `row_size` (max is {cpp.max_row_size})")
        self._n_dim = int(n_dim)
        self._row_size = int(row_size)
        self._lts = bool(local_time_stepping)
        #self._solver = new cpp.Solver(n_dim, row_size, root_mesh_size, local_time_stepping, cpp.inviscid, cpp.inviscid)
    def __dealloc__(self):
        if self._is_init:
            del self._solver

    def iteration_status(self):
        r"""! \brief fetches the `Iteration_status` describing the state of the simulation """
        status = Iteration_status()
        status.status = self._solver[0].iteration_status()
        return status

    def visualize_field(self, file_name, int n_sample = 20):
        r"""! \brief Visualizes the field data.
        \details Data is written in Tecplot subzone data format (`.szplt`).
        \param file_name name of file to write the data to, without file type extension
        \param n_sample each element will be a `n_sample`[*`n_sample`[*`n_sample]]
                        array of uniformly spaced sample points
        """
        self._solver[0].calc_jacobian()
        self._solver[0].visualize_field_tecplot(bytes(file_name, "ascii"), n_sample, False, False, True)
