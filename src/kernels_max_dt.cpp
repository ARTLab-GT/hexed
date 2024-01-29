#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed
{

#define COMPUTE_MAX_DT(Pde_templ, ...) \
{ \
  return std::min((*kernel_factory<Spatial<Pde_templ, false>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "compute time step"), \
                  (*kernel_factory<Spatial<Pde_templ,  true>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "compute time step")); \
}

double max_dt_euler(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(pde::Navier_stokes<false>::Pde)
double max_dt_navier_stokes(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time,
                            Transport_model visc, Transport_model therm_cond, bool laplacian_av)
  COMPUTE_MAX_DT(pde::Navier_stokes<true>::Pde, visc, therm_cond, laplacian_av)
double max_dt_advection(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time, double advect_length)
  COMPUTE_MAX_DT(pde::Advection, advect_length)
double max_dt_smooth_av(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(pde::Smooth_art_visc, 1., 1.)
double max_dt_fix_therm_admis(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) COMPUTE_MAX_DT(pde::Fix_therm_admis)

#undef COMPUTE_MAX_DT

}
