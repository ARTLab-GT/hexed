#include <Boundary_func.hpp>
#include <Domain_func.hpp>
#include <connection.hpp>
#include <hil_properties.hpp>

namespace hexed
{

Boundary_expr::Boundary_expr(Struct_expr expr, const Interpreter& inter) : _expr{expr}, _inter{inter} {}

std::vector<double> Boundary_expr::operator()(Boundary_connection& con, int i_fqpoint, double time) const
{
  auto sub = _inter.make_sub();
  hil_properties::surface(*sub.variables, con, i_fqpoint);
  sub.variables->assign("time", time);
  return _expr.eval(sub);
}

std::vector<double> Viscous_stress::operator()(Boundary_connection& con, int i_fqpoint, double time) const
{
  auto params = con.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  // fetch surface normal
  int nrml_sign = 1 - 2*con.inside_face_sign();
  double nrml_mag = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    double n = con.surface_normal()[i_dim*nfq + i_fqpoint];
    nrml_mag += n*n;
  }
  nrml_mag = std::sqrt(nrml_mag);
  std::vector<double> stress;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    // fetch momentum flux in reference space
    double ref_stress = -con.state_cache()[i_dim*nfq + i_fqpoint];
    // compute stress
    stress.push_back((nrml_mag > 1e-3) ? ref_stress*nrml_sign/nrml_mag : 0.);
  }
  return stress;
}

std::vector<double> Heat_flux::operator()(Boundary_connection& con, int i_fqpoint, double time) const
{
  auto params = con.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  // fetch surface normal
  int nrml_sign = 1 - 2*con.inside_face_sign();
  double nrml_mag = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    double n = con.surface_normal()[i_dim*nfq + i_fqpoint];
    nrml_mag += n*n;
  }
  nrml_mag = std::sqrt(nrml_mag);
  // fetch flux in reference space
  double ref_flux = con.ghost_face()[(params.n_dim + 1)*nfq + i_fqpoint];
  // compute flux in physical space
  return {(nrml_mag > 1e-3) ? ref_flux*nrml_sign/nrml_mag : 0.};
}

}
