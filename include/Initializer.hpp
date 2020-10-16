#ifndef CARTDG_INITIALIZER_HPP_
#define CARTDG_INITIALIZER_HPP_

#include <vector>

namespace cartdg
{

class Initializer
{
  public:
  Initializer();
  virtual ~Initializer();
  virtual std::vector<double> momentum(std::vector<double> position) = 0;
  virtual std::vector<double> scalar_state(std::vector<double> position) = 0;
};

class Constant_initializer : public Initializer
{
  public:
  Constant_initializer(int n_dim_arg, std::vector<double> state_arg);
  virtual std::vector<double> momentum(std::vector<double> position);
  virtual std::vector<double> scalar_state(std::vector<double> position);

  const int n_dim;
  std::vector<double> state;
};

}
#endif
