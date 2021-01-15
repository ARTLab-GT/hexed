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

}

#endif
