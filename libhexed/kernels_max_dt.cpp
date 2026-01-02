#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed {

#define COMPUTE_MAX_DT(Pde_templ, ...) { \
  return std::min((*kernel_factory<Spatial<Pde_templ, false>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "compute time step"), \
                  (*kernel_factory<Spatial<Pde_templ,  true>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "compute time step")); \
}

typedef pde::Navier_stokes<false, laminar> euler;
typedef pde::Navier_stokes<true, laminar> ns;
typedef pde::Navier_stokes<true, k_omega> rans;

double max_dt_euler(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(euler::Pde)
}

double max_dt_navier_stokes(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time,
                            Transport_model visc, Transport_model therm_cond) {
  bool calc_ts_ratio = true;
  if (mesh.turb_model == k_omega) COMPUTE_MAX_DT(rans::Pde, visc, therm_cond)
  else COMPUTE_MAX_DT(ns::Pde, visc, therm_cond)
}

double max_dt_advection(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time, double advect_length) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(pde::Advection, advect_length, 0)
}

double max_dt_smooth_av(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  COMPUTE_MAX_DT(pde::Smooth_art_visc, 1., 1.)
}

double max_dt_fix_nonphysical(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, bool local_time) {
  bool calc_ts_ratio = false;
  if (mesh.turb_model == k_omega) COMPUTE_MAX_DT(pde::Fix_nonphysical<4>::Pde)
  else COMPUTE_MAX_DT(pde::Fix_nonphysical<2>::Pde)
}

double max_dt_eikonal(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, double smoothing) {
  bool calc_ts_ratio = false;
  bool local_time = true;
  COMPUTE_MAX_DT(pde::Eikonal, smoothing)
  (*kernel_factory<Spatial<pde::Eikonal, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var,
                                                             smoothing))(mesh.elems);
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs(mesh.face_refinements);
  (*kernel_factory<Spatial<pde::Eikonal,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size,
                                                                  mesh.basis, mesh.mask_level, mesh.n_var))
                                                                 (face_refs, opts.sw_pr);
}

#undef COMPUTE_MAX_DT
}
