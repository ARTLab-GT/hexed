#include <Solver.hpp>
#include <Spatial.hpp>
#include <pde.hpp>
#include <Restrict_refined.hpp>

namespace hexed
{

void Solver::euler(double dt, int i_stage)
{
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  const int nd = params.n_dim;
  const int rs = params.row_size;
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes<false>::Pde>::Neighbor>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<false>::Pde>::Neighbor>(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor");
  (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
  (*kernel_factory<Spatial<Element         , pde::Navier_stokes<false>::Pde>::Local>(nd, rs, basis, dt, i_stage))(acc_mesh.cartesian().elements(), sw_car, "local");
  (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<false>::Pde>::Local>(nd, rs, basis, dt, i_stage))(acc_mesh.deformed ().elements(), sw_def, "local");
}

}
