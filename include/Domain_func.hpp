#ifndef CARTDG_DOMAIN_FUNC_HPP_
#define CARTDG_DOMAIN_FUNC_HPP_

#include "Spacetime_func.hpp"

namespace cartdg
{

class Domain_func
{
  public:
  virtual ~Domain_func();
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) = 0;
};

class State_variables : public Domain_func
{
  public:
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Error_func : public Domain_func
{
  public:
  Spacetime_func& correct;
  Error_func(Spacetime_func& correct_arg);
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

}

#endif
