#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <Prolong_refined.hpp>
#include <Restrict_refined.hpp>

namespace hexed
{

void Solver::navier_stokes(double dt, int i_stage)
{
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& bc_cons {acc_mesh.boundary_connections()};
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes>::Neighbor>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes>::Neighbor>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor");
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
  (*kernel_factory<Restrict_refined>(nd, rs, basis, false, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes>::Local>(nd, rs, basis, dt, i_stage))(acc_mesh.cartesian().elements(), sw_car, "local");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes>::Local>(nd, rs, basis, dt, i_stage))(acc_mesh.deformed ().elements(), sw_def, "local");
  (*kernel_factory<Prolong_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->apply_flux(bc_cons[i_con]);
  }
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes>::Neighbor_reconcile>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes>::Neighbor_reconcile>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor");
  (*kernel_factory<Restrict_refined>(nd, rs, basis, true, true))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage))(acc_mesh.cartesian().elements(), sw_car, "reconcile LDG flux");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes>::Reconcile_ldg_flux>(nd, rs, basis, dt, i_stage))(acc_mesh.deformed ().elements(), sw_def, "reconcile LDG flux");
}

}
