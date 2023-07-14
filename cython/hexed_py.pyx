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

cdef class Mesh:
    r"""! \brief Python interface for `hexed::Mesh` """
    cdef cpp.Mesh* mesh
    def add_tree(self, origin = np.zeros(3)):
        r"""! \see `hexed::Mesh::add_tree` """
        cdef cpp.vector[cpp.Flow_bc*] bcs
        for i_bc in range(4):
            bcs.push_back(new cpp.Nonpenetration())
        self.mesh[0].add_tree(bcs, matrix(origin))

cdef class Iteration_status:
    r"""! \brief Python interface for `hexed::Iteration_status` """
    cdef cpp.Iteration_status status
    def header(self):
        return self.status.header().decode()
    def report(self):
        return self.status.report().decode()

cdef class Solver:
    r"""! \brief Python interface for `hexed::Solver` """
    cdef cpp.Solver* _solver
    def __cinit__(self, int n_dim, int row_size, double root_mesh_size, bint local_time_stepping = False):
        self._solver = new cpp.Solver(n_dim, row_size, root_mesh_size, local_time_stepping, cpp.inviscid, cpp.inviscid)
    def __dealloc__(self):
        del self._solver
    def mesh(self):
        r"""! \brief Obtains a mutable `Mesh` object that you can interact with to generate the mesh """
        m = Mesh()
        m.mesh = &self._solver[0].mesh()
        return m
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
