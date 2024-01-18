#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>

namespace hexed
{

#define COMPUTE_CONVECTION(Pde_templ, sw) \
  auto& sw_car {(sw).children.at("cartesian")}; \
  auto& sw_def {(sw).children.at("deformed" )}; \
  const int nd = params.n_dim; \
  const int rs = params.row_size; \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(nd, rs, i_stage))(acc_mesh.cartesian().kernel_connections(), sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(nd, rs, i_stage))(acc_mesh.deformed ().kernel_connections(), sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value()))(acc_mesh.cartesian().kernel_elements(), sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(nd, rs, basis, dt, i_stage, compute_residual, allow_filter && _namespace->lookup<int>("use_filter").value()))(acc_mesh.deformed ().kernel_elements(), sw_def, "local"); \
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \

void Solver::compute_inviscid(double dt, int i_stage, bool compute_residual)
{
  bool allow_filter = true;
  COMPUTE_CONVECTION(pde::Navier_stokes<false>::Pde, stopwatch)
}

void Solver::compute_advection(double dt, int i_stage)
{
  bool allow_filter = false;
  bool compute_residual = false;
  COMPUTE_CONVECTION(pde::Advection, stopwatch.children.at("set art visc").children.at("advection"))
}

#undef COMPUTE_CONVECTION
}
