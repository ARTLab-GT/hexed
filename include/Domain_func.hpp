#ifndef CARTDG_DOMAIN_FUNC_HPP_
#define CARTDG_DOMAIN_FUNC_HPP_

#include "Qpoint_func.hpp"

namespace cartdg
{

class Spacetime_func;

/*
 * Represents a function of position, time, and flow state. That is, all the
 * variables in the mathematical problem domain. Useful for defining error functions
 * and computing integrals.
 */
class Domain_func : public Qpoint_func // can be evaluated at quadrature points, so inherits from `Qpoint_func`
{
  // the following invokes `operator()(const std::vector<double>, double, std::vector<double>)`
  // on the appropriate data at the quadrature point. Declared as private to hide the
  // technicalities of overloading inherited functions.
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint);

  public:
  virtual std::vector<double> operator()(const std::vector<double> pos, double time,
                                         const std::vector<double> state) = 0;
};

class State_variables : public Domain_func
{
  public:
  // returns `state`
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

}
#endif
