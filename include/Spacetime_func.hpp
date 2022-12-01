#ifndef HEXED_SPACETIME_FUNC_HPP_
#define HEXED_SPACETIME_FUNC_HPP_

#include <vector>
#include <array>
#include "Domain_func.hpp"

namespace hexed
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

class Position_func : public Spacetime_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim;}
  virtual inline std::string variable_name(int i_var) const {return "position" + std::to_string(i_var);}
  virtual inline std::vector<double> operator()(std::vector<double> pos, double time) const {return pos;}
};

// returns an empty vector. This is useful if you just want to visualize the position
class Empty_func : public Spacetime_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 0;}
  virtual inline std::vector<double> operator()(std::vector<double> pos, double time) const {return {};}
};

class Random_func : public Spacetime_func
{
  std::vector<double> m;
  std::vector<double> v;
  int gran;
  public:
  Random_func(std::vector<double> means, std::vector<double> variations, int granularity = 1e5);
  virtual inline int n_var(int n_dim) const {return m.size();}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

class State_from_spacetime : public Spacetime_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim + 2;}
  virtual inline std::string variable_name(int i_var) const {return "state" + std::to_string(i_var);}
};

// linear combination of components of position
class Linear : public Spacetime_func
{
  Eigen::Matrix<double, 1, Eigen::Dynamic> coefs;
  public:
  inline Linear(Eigen::VectorXd arg) : coefs{arg.transpose()} {}
  virtual inline int n_var(int n_dim) const {return 1;}
  // if size of `coefs` and `pos` don't match, truncate both to the minimum of the two sizes
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
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

// Initial consition for classic Sod problem:
// G. A. Sod. A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws. JCP 27 (1978). http://dx.doi.org/10.1016/0021-9991(78)90023-2
// Only initial condition -- not time dependent
class Sod : public State_from_spacetime
{
  public:
  double heat_rat = 1.4;
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

// multivariate Gaussian distribution in terms of spatial coordinates
// normalized to evaluate to 1 at the zero vector
class Spatial_gaussian : public Spacetime_func
{
  std::vector<double> dev;
  public:
  // specify the standard deviation for each variable
  // trailing variables can be left unspecified and default to the last specified component of std_dev
  Spatial_gaussian(std::vector<double> std_dev);
  inline int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

// steady-state solution to Laplacian diffusion in an anular domain (for verification testing)
// velocity is uniformly zero, total energy is uniform, and mass is proportional to log of distance from origin
class Annular_diffusion_test : public State_from_spacetime
{
  double ref_mass;
  double ref_radius;
  double ener;
  public:
  // `reference_mass` is the value of mass at `reference_radius`
  Annular_diffusion_test(double reference_mass, double reference_radius, double energy);
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const;
};

}
#endif
