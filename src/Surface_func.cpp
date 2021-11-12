#include <cmath>

#include <Surface_func.hpp>

namespace cartdg
{

Surface_from_domain::Surface_from_domain(Domain_func& df) : domain{df} {}

std::vector<double> Surface_from_domain::operator()(std::vector<double> pos, double time,
                                                    std::vector<double> state, std::vector<double> normal)
{
  return domain(pos, time, state);
}

std::vector<double> Force_per_area::operator()(std::vector<double> pos, double time,
                                               std::vector<double> state, std::vector<double> normal)
{
  double momentum_sq {0.};
  double normal_sq {0.};
  for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
    momentum_sq += state[i_dim]*state[i_dim];
    normal_sq += normal[i_dim]*normal[i_dim];
  }
  double normal_mag {std::sqrt(normal_sq)};
  double pres {0.4*(state[state.size() - 1] - 0.5*momentum_sq/state[state.size() - 2])};
  std::vector<double> fpa;
  for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
    fpa.push_back(-normal[i_dim]*pres/normal_mag);
  }
  return fpa;
}

}
