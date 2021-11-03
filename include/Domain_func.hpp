#ifndef CARTDG_DOMAIN_FUNC_HPP_
#define CARTDG_DOMAIN_FUNC_HPP_

#include "Spacetime_func.hpp"

namespace cartdg
{

class Domain_func
{
  public:
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) = 0;
};

class State_variables : public Domain_func
{
  public:
  // returns `state`
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Domain_func_from_st : public Domain_func
{
  Spacetime_func& spacetime;
  public:
  Domain_func_from_st(Spacetime_func&);
  // evaluates given Spacetime_func at `(point_pos, point_time)`
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Error_func : public Domain_func
{
  Spacetime_func& correct;
  public:
  Error_func(Spacetime_func&);
  // returns elementwise difference between `state` and `correct(point_pos, point_time)`
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

}

#endif
