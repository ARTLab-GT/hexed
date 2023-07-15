from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Eigen/Dense" namespace "Eigen":
    cdef cppclass MatrixXd:
        MatrixXd()
        MatrixXd(int rows, int cols)
        int size()
        int rows()
        int cols()
        double* data()
        void resize(int rows, int cols)

cdef extern from "config.hpp" namespace "hexed::config":
    cdef int max_row_size

cdef extern from "Boundary_condition.hpp" namespace "hexed":
    cdef cppclass Flow_bc:
        pass
    cdef cppclass Nonpenetration(Flow_bc):
        pass

cdef extern from "Mesh.hpp" namespace "hexed":
    cdef cppclass Mesh:
        void add_tree(vector[Flow_bc*] extremal_bcs, MatrixXd origin) except+

cdef extern from "Iteration_status.hpp" namespace "hexed":
    cdef cppclass Iteration_status:
        string header() except+
        string report() except+

cdef extern from "Transport_model.hpp" namespace "hexed":
    cdef cppclass Transport_model:
        const bool is_viscous
        double coefficient(double sqrt_temp) const
        @staticmethod
        Transport_model inviscid() except+
        @staticmethod
        Transport_model constant(double value) except+
        @staticmethod
        Transport_model sutherland(double reference_value, double reference_temperature, double temperature_offset) except+
    const Transport_model inviscid
    const Transport_model air_const_dyn_visc
    const Transport_model air_const_therm_cond
    const Transport_model air_sutherland_dyn_visc
    const Transport_model air_sutherland_therm_cond

cdef extern from "Solver.hpp" namespace "hexed":
    cdef cppclass Solver:
        Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping, Transport_model viscosity_model, Transport_model thermal_conductivity_model) except+
        Iteration_status iteration_status() except+
        Mesh& mesh() except+
        void calc_jacobian() except+
        void visualize_field_tecplot(string name, int n_sample, bool edges, bool qpoints, bool interior) except+
