#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <Prolong_refined.hpp>
#include <Restrict_refined.hpp>

namespace hexed
{

#define COMPUTE_DIFFUSION(Pde_templ, sw, bc_fun) \
  auto& sw_car {(sw).children.at("cartesian")}; \
  auto& sw_def {(sw).children.at("deformed" )}; \
  const int nd = params.n_dim; \
  const int rs = params.row_size; \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(nd, rs, i_stage PDE_ARGS))(acc_mesh.cartesian().kernel_connections(), sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(nd, rs, i_stage PDE_ARGS))(acc_mesh.deformed ().kernel_connections(), sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis, false, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value() PDE_ARGS))(acc_mesh.cartesian().kernel_elements(), sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value() PDE_ARGS))(acc_mesh.deformed ().kernel_elements(), sw_def, "local"); \
  if (!i_stage) { \
    (*kernel_factory<Prolong_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
    bc_fun(); \
    (*kernel_factory<Spatial<Pde_templ, false>::Neighbor_reconcile>(nd, rs))(acc_mesh.cartesian().kernel_connections(), sw_car, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor_reconcile>(nd, rs))(acc_mesh.deformed ().kernel_connections(), sw_def, "neighbor"); \
    (*kernel_factory<Restrict_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
    (*kernel_factory<Spatial<Pde_templ, false>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value()))(acc_mesh.cartesian().kernel_elements(), sw_car, "reconcile LDG flux"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value()))(acc_mesh.deformed ().kernel_elements(), sw_def, "reconcile LDG flux"); \
  } \
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \

#define PDE_ARGS , visc, therm_cond, bool(_namespace->lookup<int>("laplacian_art_visc").value())
void Solver::compute_viscous(double dt, int i_stage, bool compute_residual)
{
  bool allow_filter = true;
  COMPUTE_DIFFUSION(pde::Navier_stokes<true>::Pde, stopwatch, apply_flux_bcs)
}
#undef PDE_ARGS

#define PDE_ARGS
void Solver::compute_fta(double dt, int i_stage)
{
  bool allow_filter = false;
  bool compute_residual = false;
  COMPUTE_DIFFUSION(pde::Fix_therm_admis, stopwatch.children.at("fix admis."), apply_fta_flux_bcs)
}

void Solver::compute_avc_diff(double dt, int i_stage)
{
  bool allow_filter = false;
  bool compute_residual = false;
  COMPUTE_DIFFUSION(pde::Smooth_art_visc, stopwatch.children.at("set art visc").children.at("diffusion"), apply_avc_diff_flux_bcs)
}
#undef PDE_ARGS

#undef COMPUTE_DIFFUSION
}
