from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr

cdef extern from "Eigen/Dense" namespace "Eigen":
    cdef cppclass MatrixXd:
        MatrixXd()
        MatrixXd(int rows, int cols)
        int size()
        int rows()
        int cols()
        double* data()
        void resize(int rows, int cols)
        double& operator()(int row, int col)
    cdef cppclass Matrix2d:
        int size()
        int rows()
        int cols()
        double* data()
        double& operator()(int row, int col)
    cdef cppclass MatrixXi:
        MatrixXi(int rows, int cols)
        int size()
        int rows()
        int cols()
        int* data()
        int& operator()(int row, int col)

cdef extern from "config.hpp" namespace "hexed::config":
    cdef int max_row_size

cdef extern from "Boundary_condition.hpp" namespace "hexed":
    cdef cppclass Flow_bc:
        pass
    cdef cppclass Nonpenetration(Flow_bc):
        pass

cdef extern from "Surface_geom.hpp" namespace "hexed":
    cdef cppclass Surface_geom:
        pass
cdef extern from "Simplex_geom.hpp" namespace "hexed":
    cdef cppclass Simplex_geom2(Surface_geom):
        Simplex_geom2(vector[Matrix2d]&&) except+
    cdef cppclass Simplex_geom3(Surface_geom):
        Simplex_geom3(vector[MatrixXd]&&) except+
    vector[Matrix2d] segments(const MatrixXd& points) except+
cdef extern from "Element.hpp" namespace "hexed":
    cdef cppclass Element:
        double resolution_badness

ctypedef bool(*ref_criterion)(Element&)

cdef extern from "Mesh.hpp" namespace "hexed::Mesh":
    cdef cppclass General_ref_criterion:
        General_ref_criterion(double, int, int)
cdef extern from "Mesh.hpp" namespace "hexed":
    cdef cppclass Mesh:
        void add_tree(vector[Flow_bc*] extremal_bcs, MatrixXd origin) except+
        void update() except+
        void update(General_ref_criterion) except+
        void set_surface(Surface_geom*, Flow_bc*, MatrixXd) except+
        void relax() except+
        int surface_bc_sn()
        @staticmethod
        bool always(Element&)
        @staticmethod
        bool never(Element&)
        @staticmethod
        bool if_extruded(Element&)
        @staticmethod
        ref_criterion Res_bad_tol(double tol)

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
        void set_res_bad_surface_rep(int bc_sn) except+
        void calc_jacobian() except+
        void visualize_field_tecplot(string name, int n_sample, bool edges, bool qpoints, bool interior) except+
