#ifndef CARTDG_DOMAIN_FUNC_HPP_
#define CARTDG_DOMAIN_FUNC_HPP_

#include "Spacetime_func.hpp"

namespace cartdg
{

class Domain_func
{
  public:
  virtual ~Domain_func() = default;
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

class Domain_from_spacetime : public Domain_func
{
  Spacetime_func& spacetime;
  public:
  Domain_from_spacetime(Spacetime_func&);
  // evaluates given Spacetime_func at `(point_pos, point_time)`
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Diff_sq : public Domain_func
{
  Domain_func& func0;
  Domain_func& func1;
  public:
  // arguments must return values of same size
  Diff_sq(Domain_func&, Domain_func&);
  // returns elementwise squared difference between provided funcs
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Error_func : public Domain_func
{
  Spacetime_func& correct;
  Domain_from_spacetime dfs;
  State_variables sv;
  Diff_sq ds;
  public:
  Error_func(Spacetime_func&);
  // returns elementwise difference between `state` and `correct(point_pos, point_time)`, squared
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Stag_pres : public Domain_func
{
  double hr;
  public:
  Stag_pres(double heat_rat = 1.4);
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

class Stag_pres_errsq : public Domain_func
{
  Stag_pres sp;
  Constant_func free;
  Domain_from_spacetime dfs;
  Diff_sq ds;
  public:
  Stag_pres_errsq(std::vector<double> freestream, double heat_rat = 1.4);
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state);
};

}

#endif
