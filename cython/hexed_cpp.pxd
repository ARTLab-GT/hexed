from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "Iteration_status.hpp" namespace "hexed":
    cdef cppclass Iteration_status:
        string header()
        string report()

cdef extern from "Transport_model.hpp" namespace "hexed":
    cdef cppclass Transport_model:
        const bool is_viscous
        double coefficient(double sqrt_temp) const
        @staticmethod
        Transport_model inviscid()
        @staticmethod
        Transport_model constant(double value)
        @staticmethod
        Transport_model sutherland(double reference_value, double reference_temperature, double temperature_offset)
    const Transport_model inviscid
    const Transport_model air_const_dyn_visc
    const Transport_model air_const_therm_cond
    const Transport_model air_sutherland_dyn_visc
    const Transport_model air_sutherland_therm_cond

cdef extern from "Solver.hpp" namespace "hexed":
    cdef cppclass Solver:
        Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping, Transport_model viscosity_model, Transport_model thermal_conductivity_model)
        Iteration_status iteration_status()
