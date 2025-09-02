#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed {

#define COMPUTE_MAX_DT(Pde_templ, ...) { \
  return std::min((*kernel_factory<Spatial<Pde_templ, false>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "compute time step"), \
                  (*kernel_factory<Spatial<Pde_templ,  true>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "compute time step")); \
}

typedef pde::Navier_stokes<false, laminar> euler;
typedef pde::Navier_stokes<true, laminar> ns;
typedef pde::Navier_stokes<true, k_omega> rans;
double max_dt_euler(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(euler::Pde)
double max_dt_navier_stokes(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time,
                            Transport_model visc, Transport_model therm_cond) {
  if (mesh.turb_model == k_omega) COMPUTE_MAX_DT(rans::Pde, visc, therm_cond)
  else COMPUTE_MAX_DT(ns::Pde, visc, therm_cond)
}
double max_dt_advection(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time, double advect_length)
  COMPUTE_MAX_DT(pde::Advection, advect_length)
double max_dt_smooth_av(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(pde::Smooth_art_visc, 1., 1.)
double max_dt_fix_therm_admis(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(pde::Fix_therm_admis)

#undef COMPUTE_MAX_DT

}
