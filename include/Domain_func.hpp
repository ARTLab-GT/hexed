#ifndef CARTDG_DOMAIN_FUNC_HPP_
#define CARTDG_DOMAIN_FUNC_HPP_

#include <vector>

#include "Spacetime_func.hpp"

namespace cartdg
{

class Domain_func
{
  public:
  virtual double operator()(const std::vector<double> point_pos, double point_time,
                            const std::vector<double> state) = 0;
};

class State_variable : public Domain_func
{
  public:
  int i_var;
  State_variable();
  State_variable(int i_var_arg);
  virtual double operator()(const std::vector<double> point_pos, double point_time,
                            const std::vector<double> state);
};

class Error_func : public Domain_func
{
  public:
  std::vector<double> weights;
  Spacetime_func& correct;
  Error_func(Spacetime_func& correct_arg);
  Error_func(Spacetime_func& correct_arg, std::vector<double> weights_arg);
  virtual double operator()(const std::vector<double> point_pos, double point_time,
                            const std::vector<double> state);
};

}

#endif
