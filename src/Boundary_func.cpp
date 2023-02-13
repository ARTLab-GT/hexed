#include <Boundary_func.hpp>
#include <Domain_func.hpp>

namespace hexed
{

std::vector<double> Viscous_stress::operator()(Boundary_face& bf, int i_fqpoint, double time) const
{
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  std::vector<double> stress;
  int nrml_sign = 1 - 2*bf.inside_face_sign();
  double nrml_mag = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    double n = bf.surface_normal()[i_dim*nfq + i_fqpoint];
    nrml_mag += n*n;
  }
  nrml_mag = std::sqrt(nrml_mag);
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    stress.push_back(-bf.inside_face()[(2*params.n_var + i_dim)*nfq + i_fqpoint]*nrml_sign/nrml_mag);
  }
  return stress;
}

std::vector<double> Heat_flux::operator()(Boundary_face& bf, int i_fqpoint, double time) const
{
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  int nrml_sign = 1 - 2*bf.inside_face_sign();
  double nrml_mag = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    double n = bf.surface_normal()[i_dim*nfq + i_fqpoint];
    nrml_mag += n*n;
  }
  nrml_mag = std::sqrt(nrml_mag);
  return {bf.inside_face()[(2*params.n_var + params.n_dim + 1)*nfq + i_fqpoint]*nrml_sign/nrml_mag};
}

}
