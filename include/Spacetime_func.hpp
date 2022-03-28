#ifndef CARTDG_SPACETIME_FUNC_HPP_
#define CARTDG_SPACETIME_FUNC_HPP_

#include <vector>
#include <array>
#include "Domain_func.hpp"

namespace cartdg
{

/*
 * Represents a function of position and time. Useful for specifying initial conditions
 * and analytic solutions.
 */
class Spacetime_func : public Domain_func
{
  // The following invokes `operator()(pos, time)`. Declared as private to avoid the
  // technicalities of overloading inherited functions.
  virtual std::vector<double> operator()(const std::vector<double> pos, double time,
                                         const std::vector<double> state) const;
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const = 0;
};

class Constant_func : public Spacetime_func
{
  std::vector<double> value;
  public:
  Constant_func(std::vector<double> value_arg);
  virtual inline int n_var(int n_dim) const {return value.size();}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const; // returns `value`
};

class State_from_spacetime : public Spacetime_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim + 2;}
  virtual inline std::string variable_name(int i_var) const {return "state" + std::to_string(i_var);}
};

/*
 * Isentropic vortex flow described by Gaussian function, a classic
 * CFD test problem. Flow field is an exact solution of the Euler equations
 * which consists of the vortex translating at the freestream velocity. Flowfield
 * is well defined at all position and all time.
 */
class Isentropic_vortex : public State_from_spacetime
{
  std::vector<double> freestream;

  public:
  double heat_rat = 1.4;
  double argmax_radius = 0.05; // velocity reaches its maximum value at this radius
  double max_nondim_veloc = 0.02; // max velocity perturbation normalized by freestream speed of sound
  double center0 = 0.;
  double center1 = 0.;
  Isentropic_vortex(std::vector<double> freestream_state);
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

/*
 * Velocity field of this flow matches an irrotational, incompressible source/vortex doublet.
 * The circle of radius `radius` will be a streamline of this flow. Thermodynamic
 * variables are set such that entropy and stagnation enthalpy are constant.
 * flow is isentropic. Flow field is steady, and is an exact solution only in
 * incompressible flow.
 *
 * WARNING: Singularity at `location`! This class is applicable only to domains
 * which do not include this point.
 */
class Doublet : public State_from_spacetime
{
  std::vector<double> freestream;
  int n_v;
  int n_dim;

  public:
  std::array<double, 2> location {0., 0.};
  double radius {1.};
  double heat_rat {1.4};
  Doublet(std::vector<double> freestream_state);
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

}
#endif
