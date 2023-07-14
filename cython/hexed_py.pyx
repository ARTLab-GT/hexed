import numpy as np
cimport hexed_cpp as cpp

## \namespace hexed_py
# \brief foo

cdef matrix_shape(arr):
    """! @brief gets the shape of a matrix """
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
    cdef cpp.Mesh* mesh
    def add_tree(self, origin = np.zeros(3)):
        cdef cpp.vector[cpp.Flow_bc*] bcs
        for i_bc in range(4):
            bcs.push_back(new cpp.Nonpenetration())
        self.mesh[0].add_tree(bcs, matrix(origin))

cdef class Iteration_status:
    cdef cpp.Iteration_status status
    def header(self):
        return self.status.header().decode()
    def report(self):
        return self.status.report().decode()

cdef class Solver:
    cdef cpp.Solver* _solver
    def __cinit__(self, int n_dim, int row_size, double root_mesh_size, bint local_time_stepping = False):
        self._solver = new cpp.Solver(n_dim, row_size, root_mesh_size, local_time_stepping, cpp.inviscid, cpp.inviscid)
    def __dealloc__(self):
        del self._solver
    def mesh(self):
        m = Mesh()
        m.mesh = &self.sol[0].mesh()
        return m
    def iteration_status(self):
        status = Iteration_status()
        status.status = self._solver[0].iteration_status()
        return status
    def visualize(self, file_name, int n_sample = 20):
        self._solver[0].calc_jacobian()
        self._solver[0].visualize_field_tecplot(bytes(file_name, "ascii"), n_sample, False, False, True)
