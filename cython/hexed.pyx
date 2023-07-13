cimport hexed_cpp as cpp

cdef class Iteration_status:
    cdef cpp.Iteration_status status
    @staticmethod
    cdef make(cpp.Iteration_status arg):
        made = Iteration_status()
        made.status = arg
        return made
    def header(self):
        return str(self.status.header())
    def report(self):
        return str(self.status.report())

cdef class Solver:
    cdef cpp.Solver* sol
    def __cinit__(self, int n_dim, int row_size, double root_mesh_size, bint local_time_stepping = False):
        self.sol = new cpp.Solver(n_dim, row_size, root_mesh_size, local_time_stepping, cpp.inviscid, cpp.inviscid)
    def __dealloc__(self):
        del self.sol
    def iteration_status(self):
        return Iteration_status.make(self.sol[0].iteration_status())
