#ifndef HEXED_SPACETIME_FUNC_HPP_
#define HEXED_SPACETIME_FUNC_HPP_

#include <vector>
#include <array>
#include "Domain_func.hpp"
#include "config.hpp"

namespace hexed
{

/*!
 * Represents a function of position and time.
 * Useful for specifying initial conditions and analytic flow solutions.
 */
class Spacetime_func : public Domain_func
{
  std::vector<double> operator()(std::vector<double> pos, double time,
                                 std::vector<double> state) const override;
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const = 0;
};

/*! \brief Evaluates a \ref struct_expr "structured expression"
 * \details Expression is evaluated in an environment that includes
 * `pos0`, `pos1`, `pos2`, and `time`.
 */
class Spacetime_expr : public Spacetime_func
{
  Struct_expr _expr;
  const Interpreter& _inter;
  public:
  Spacetime_expr(Struct_expr, const Interpreter&);
  Spacetime_expr(Struct_expr, Interpreter&&) = delete;
  inline int n_var(int n_dim) const override {return _expr.names.size();}
  inline std::string variable_name(int n_dim, int i_var) const override {return _expr.names[i_var];}
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

//! Always returns the same constant (vector) value.
class Constant_func : public Spacetime_func
{
  std::vector<double> value;
  public:
  //! \param value_arg The constant value to be returned. Size can be whatever you want.
  Constant_func(std::vector<double> value_arg);
  inline int n_var(int n_dim) const override {return value.size();}
  inline std::string variable_name(int n_dim, int i_var) const override {return "constant" + std::to_string(i_var);}
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

//! Returns the position vector.
class Position_func : public Spacetime_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "position" + std::to_string(i_var);}
  inline std::vector<double> operator()(std::vector<double> pos, double time) const override {return pos;}
};

//! Returns an empty vector. You can pass this to `Solver::visualize_field_tecplot` if you just want to visualize the position.
class Empty_func : public Spacetime_func
{
  public:
  inline int n_var(int n_dim) const override {return 0;}
  inline std::vector<double> operator()(std::vector<double> pos, double time) const override {return {};}
};

//! Returns a random output uniformly distributed in a user-specified range.
class Random_func : public Spacetime_func
{
  std::vector<double> m;
  std::vector<double> v;
  int gran;
  public:
  /*!
   * \param means The mean of each variable. The size of `means` determines the number of output variables.
   * \param variations The size of the range for each variable. Must be the same size as `means`.
   * \param granularity Results will be divisible by `1./granularity`.
   * \details Output variable `i` will be in the range `[means[i] - variations[i]/2, means[i] + variations[i]/2]`.
   */
  Random_func(std::vector<double> means, std::vector<double> variations, int granularity = 1e5);
  inline int n_var(int n_dim) const override {return m.size();}
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

//! Computes a single output variable which is a linear combination of components of position.
class Linear : public Spacetime_func
{
  Eigen::Matrix<double, 1, Eigen::Dynamic> coefs;
  public:
  /*! \param arg Coefficients for components of position.
   * \details If `arg.size() > n_dim`, only the first `n_dim` components of `arg` will be used.
   * If `arg.size() < n_dim`, only the first `arg.size()` components of position will be used.
   */
  inline Linear(Eigen::VectorXd arg) : coefs{arg.transpose()} {}
  inline int n_var(int n_dim) const override {return 1;}
  // if size of `coefs` and `pos` don't match, truncate both to the minimum of the two sizes
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

//! A class of `Spacetime_func`s whose output is a state vector.
class State_from_spacetime : public Spacetime_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim + 2;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "state" + std::to_string(i_var);}
};

/*!
 * Velocity field of this flow matches an irrotational, incompressible source/vortex doublet.
 * The circle of radius `radius` will be a streamline of this flow. Thermodynamic
 * variables are set such that entropy and stagnation enthalpy are constant.
 * flow is isentropic. Flow field is steady, and is an exact solution only in
 * incompressible flow.
 *
 * \attention Singularity at `location`! This class is applicable only to domains
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
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

/*! \brief Initial consition for classic Sod problem.
 * \details G. A. Sod. A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws. JCP 27 (1978). http://dx.doi.org/10.1016/0021-9991(78)90023-2
 * \attention Only initial condition -- not time dependent!
 */
class Sod : public State_from_spacetime
{
  public:
  double heat_rat = 1.4;
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

//! multivariate Gaussian distribution in terms of spatial coordinates normalized to evaluate to 1 at the zero vector
class Spatial_gaussian : public Spacetime_func
{
  std::vector<double> dev;
  public:
  /*! \param std_dev the standard deviation for each variable.
   * Trailing variables can be left unspecified and default to the last specified component of `std_dev`.
   */
  Spatial_gaussian(std::vector<double> std_dev);
  inline int n_var(int n_dim) const {return 1;}
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

/*! \brief Steady-state solution to Laplacian diffusion in an anular domain (for verification testing).
 * \details Velocity is uniformly zero, total energy is uniform, and mass is proportional to log of distance from origin.
 */
class Annular_diffusion_test : public State_from_spacetime
{
  double val_scale;
  double rad_scale;
  double ener;
  public:
  //! \verbatim mass = value_scalar*std::log(radius/radius_scalar) \endverbatim
  Annular_diffusion_test(double value_scalar, double radius_scalar, double energy);
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};

#if HEXED_USE_NLOPT
/*! \brief Steady-state solution for Ringleb flow.
 * \details [Ringleb flow](https://www.cfd-online.com/Wiki/Ringleb_flow)
 * is a steady, sub/transonic, isentropic inviscid flow.
 * Exact solution is available as a closed form for position
 * in terms of velocity and stream function.
 * This `Spacetime_func` gives you state as a function of position using a numerical root finder,
 * so be sure to set the tolerance if you care where you stand on the accuracy/speed tradeoff.
 */
class Ringleb : public State_from_spacetime
{
  double tol;
  double heat_rat;
  public:
  inline Ringleb(double root_tolerance = 1e-12, double heat_ratio = 1.4) : tol{root_tolerance}, heat_rat{heat_ratio} {}
  std::vector<double> operator()(std::vector<double> pos, double time) const override;
};
#endif

}
#endif
