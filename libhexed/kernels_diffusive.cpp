#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed {

#define COMPUTE_DIFFUSION(Pde_templ, ...) { \
  Vector_view<Hard_kernel_connection> car_cons(mesh.car_connections); \
  Vector_view<Hard_kernel_connection> def_cons(mesh.def_connections); \
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs(mesh.face_refinements); \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, false, true))(face_refs, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  if (!opts.i_stage) { \
    (*kernel_factory<Spatial<Pde_templ, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, true, true))(face_refs, opts.sw_pr); \
    flux_bc(); \
    (*kernel_factory<Spatial<Pde_templ, false>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size, mesh.mask_level))(car_cons, opts.sw_car, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size, mesh.mask_level))(def_cons, opts.sw_def, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, true, true))(face_refs, opts.sw_pr); \
    (*kernel_factory<Spatial<Pde_templ, false>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "reconcile LDG flux"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "reconcile LDG flux"); \
  } \
  (*kernel_factory<Spatial<Pde_templ,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
}

typedef pde::Navier_stokes<true, laminar> ns;
typedef pde::Navier_stokes<true, k_omega> rans;
void compute_navier_stokes(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc, Transport_model visc, Transport_model therm_cond, bool update_prod) {
  if (mesh.turb_model == k_omega) COMPUTE_DIFFUSION(rans::Pde, visc, therm_cond)
  else COMPUTE_DIFFUSION(ns::Pde, visc, therm_cond)
}
void compute_smooth_av(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc, double diff_time, double cheby_step) {
  COMPUTE_DIFFUSION(pde::Smooth_art_visc, diff_time, cheby_step)
}
void compute_fix_therm_admis(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc) {
  COMPUTE_DIFFUSION(pde::Fix_therm_admis)
}

#undef COMPUTE_DIFFUSION

}
