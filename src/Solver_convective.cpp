#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>

namespace hexed
{

#define LOCAL(Pde_templ, sw, axi) { \
  (*kernel_factory<Spatial<Element         , Pde_templ, axi>::Local>(nd, rs, basis, dt, i_stage, _namespace->lookup<int>("use_filter").value()))(acc_mesh.cartesian().elements(), sw_car, "local"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ, axi>::Local>(nd, rs, basis, dt, i_stage, _namespace->lookup<int>("use_filter").value()))(acc_mesh.deformed ().elements(), sw_def, "local"); \
}

#define COMPUTE_CONVECTION(Pde_templ, sw) \
  auto& sw_car {(sw).children.at("cartesian")}; \
  auto& sw_def {(sw).children.at("deformed" )}; \
  const int nd = params.n_dim; \
  const int rs = params.row_size; \
  (*kernel_factory<Spatial<Element         , Pde_templ>::Neighbor>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Deformed_element, Pde_templ>::Neighbor>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor"); \
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \
  if (_namespace->lookup<int>("axisymmetric").value()) LOCAL(Pde_templ, sw, true) \
  else LOCAL(Pde_templ, sw, false) \
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict")); \

void Solver::compute_inviscid(double dt, int i_stage)
{
  COMPUTE_CONVECTION(pde::Navier_stokes<false>::Pde, stopwatch)
}

void Solver::compute_advection(double dt, int i_stage)
{
  COMPUTE_CONVECTION(pde::Advection, stopwatch.children.at("set art visc").children.at("advection"))
}

#undef COMPUTE_CONVECTION
#undef LOCAL
}
