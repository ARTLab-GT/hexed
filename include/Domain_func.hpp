#ifndef HEXED_DOMAIN_FUNC_HPP_
#define HEXED_DOMAIN_FUNC_HPP_

#include "Qpoint_func.hpp"
#include "Surface_func.hpp"

namespace hexed
{

class Spacetime_func;

/*
 * Represents a function of position, time, and flow state. That is, all the
 * variables in the mathematical problem domain. Useful for defining error functions
 * and computing integrals.
 */
class Domain_func : public Qpoint_func, public Surface_func
{
  // the following invoke `operator()(const std::vector<double>, double, std::vector<double>)`
  // on the appropriate data at the quadrature point. Declared as private to hide the
  // technicalities of overloading inherited functions.
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const;

  public:
  virtual std::vector<double> operator()(const std::vector<double> pos, double time,
                                         const std::vector<double> state) const = 0;
};

class State_variables : public Domain_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim + 2;}
  virtual inline std::string variable_name(int i_var) const {return "state" + std::to_string(i_var);}
  // returns `state`
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Diff_sq : public Domain_func
{
  const Domain_func& func0;
  const Domain_func& func1;
  public:
  // arguments must return values of same size
  Diff_sq(const Domain_func& , const Domain_func& ); // note: don't use const references because they would allow invalid passing of temporary object
  Diff_sq(      Domain_func&&, const Domain_func& ) = delete; // can't accept temporaries or else dangling reference
  Diff_sq(const Domain_func& ,       Domain_func&&) = delete;
  Diff_sq(      Domain_func&&,       Domain_func&&) = delete;
  virtual inline int n_var(int n_dim) const {return func0.n_var(n_dim);}
  // returns elementwise squared difference between provided funcs
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Error_func : public Domain_func
{
  const Spacetime_func& correct;
  public:
  Error_func(const Spacetime_func&);
  Error_func(Spacetime_func&&) = delete;
  virtual int n_var(int n_dim) const;
  virtual std::string variable_name(int i_var) const;
  // returns elementwise difference between `state` and `correct(point_pos, point_time)`, squared
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Stag_pres : public Domain_func
{
  double hr;
  public:
  Stag_pres(double heat_rat = 1.4);
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "stagnation_pressure";}
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Pressure : public Domain_func
{
  double hr;
  public:
  Pressure(double heat_rat = 1.4);
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "pressure";}
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Velocity : public Domain_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual std::string variable_name(int i_var) const;
  virtual std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                         const std::vector<double> state) const;
};

class Mass : public Domain_func
{
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "mass";}
  virtual inline std::vector<double> operator()(const std::vector<double> point_pos, double point_time,
                                                const std::vector<double> state) const
  {
    return {state[point_pos.size()]};
  }
};

}
#endif
