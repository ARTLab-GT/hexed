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
  (*kernel_factory<Spatial<Element         , Pde_templ>::Neighbor>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ>::Neighbor>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis, false, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Spatial<Element         , Pde_templ>::Local>(nd, rs, basis, dt, i_stage PDE_ARGS))(acc_mesh.cartesian().elements(), sw_car, "local"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ>::Local>(nd, rs, basis, dt, i_stage PDE_ARGS))(acc_mesh.deformed ().elements(), sw_def, "local"); \
  (*kernel_factory<Prolong_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  bc_fun(); \
  (*kernel_factory<Spatial<Element         , Pde_templ>::Neighbor_reconcile>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ>::Neighbor_reconcile>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Spatial<Element         , Pde_templ>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage))(acc_mesh.cartesian().elements(), sw_car, "reconcile LDG flux"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage))(acc_mesh.deformed ().elements(), sw_def, "reconcile LDG flux"); \
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \

#define PDE_ARGS , visc
void Solver::compute_viscous(double dt, int i_stage)
{
  COMPUTE_DIFFUSION(pde::Navier_stokes<true>::Pde, stopwatch, apply_flux_bcs)
}
#undef PDE_ARGS

#define PDE_ARGS
void Solver::compute_fta(double dt, int i_stage)
{
  COMPUTE_DIFFUSION(pde::Fix_therm_admis, stopwatch.children.at("fix admis."), apply_fta_flux_bcs)
}

void Solver::compute_avc_diff(double dt, int i_stage)
{
  COMPUTE_DIFFUSION(pde::Smooth_art_visc, stopwatch.children.at("set art visc").children.at("diffusion"), apply_avc_diff_flux_bcs)
}
#undef PDE_ARGS

#undef COMPUTE_DIFFUSION
}
