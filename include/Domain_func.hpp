#ifndef CARTDG_DOMAIN_FUNCTION_HPP_
#define CARTDG_DOMAIN_FUNCTION_HPP_

#include <vector>

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

}

#endif
