#include <cmath>

#include <Surface_func.hpp>
#include <Domain_func.hpp>
#include <Boundary_condition.hpp>

namespace hexed
{

std::vector<double> Surface_func::operator()(Boundary_face& bf, int i_fqpoint, double time) const
{
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  std::vector<double> pos;
  std::vector<double> normal;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    pos.push_back(bf.surface_position()[i_dim*nfq + i_fqpoint]);
    normal.push_back(bf.surface_normal()[i_dim*nfq + i_fqpoint]);
  }
  std::vector<double> state;
  for (int i_var = 0; i_var < params.n_var; ++i_var) state.push_back(bf.inside_face()[i_var*nfq + i_fqpoint]);
  return operator()(pos, time, state, normal);
}

Force_per_area::Force_per_area(double heat_rat) : hr{heat_rat} {}

std::vector<double> Force_per_area::operator()(std::vector<double> pos, double time,
                                               std::vector<double> state, std::vector<double> normal) const
{
  double momentum_sq {0.};
  double normal_sq {0.};
  for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
    momentum_sq += state[i_dim]*state[i_dim];
    normal_sq += normal[i_dim]*normal[i_dim];
  }
  double normal_mag {std::sqrt(normal_sq)};
  double pres {(hr - 1.)*(state[state.size() - 1] - 0.5*momentum_sq/state[state.size() - 2])}; // FIXME: use correct heat_rat
  std::vector<double> fpa;
  for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
    fpa.push_back(-normal[i_dim]*pres/normal_mag);
  }
  return fpa;
}

}
