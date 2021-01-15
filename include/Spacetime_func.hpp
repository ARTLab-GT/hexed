#ifndef SPACETIME_FUNC_HPP_
#define SPACETIME_FUNC_HPP_

#include <vector>

namespace cartdg
{

class Spacetime_func
{
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time) = 0;
};

class Constant_func : public Spacetime_func
{
  public:
  std::vector<double> value;
  virtual std::vector<double> operator()(std::vector<double> pos, double time);
};

class Isentropic_vortex : public Spacetime_func
{
  public:
  std::vector<double> state;
  double heat_rat = 1.4;
  double critical_radius = 0.05;
  Isentropic_vortex(std::vector<double> state_arg);
  virtual std::vector<double> operator()(std::vector<double> pos, double time);
};

}

#endif
