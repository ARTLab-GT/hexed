#include <Initializer.hpp>

namespace cartdg
{

Initializer::Initializer() {}
Initializer::~Initializer() {}

Constant_initializer::Constant_initializer(int n_dim_arg, std::vector<double> state_arg)
: n_dim(n_dim_arg), state(state_arg) {}

std::vector<double> Constant_initializer::momentum(std::vector<double> position)
{
  std::vector<double> value;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    value.push_back(state[i_dim]);
  }
  return value;
}

std::vector<double> Constant_initializer::scalar_state(std::vector<double> position)
{
  std::vector<double> value;
  for (int i_var = n_dim; i_var < (int)state.size(); ++i_var)
  {
    value.push_back(state[i_var]);
  }
  return value;
}

}
