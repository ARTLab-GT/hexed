#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include <vector>
#include "Output_data.hpp"
#include "Element.hpp"
#include "Basis.hpp"

namespace cartdg
{

class Grid;

/*
 * Represents a function which can be evaluated at quadrature points. Can depend on
 * flow state, position, time, or on mathematical parameters like element Jabobian,
 * quadrature weights, etc. Note: `Output_data` is a virtual base class because
 * derived classes may inherit from multiple types of `Output_data`.
 */
class Qpoint_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const = 0;
};

// Returns a vector with one element: the Jacobian determinant at the quadrature point.
class Jacobian_det_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "jacobian_determinant";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

class Time_step_scale_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "time_step_scale";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

class Viscosity_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "artificial_visc_coef";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

// returns the most recent update to the state divided by the local time step scale.
class Physical_update : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim + 2;}
  virtual inline std::string variable_name(int i_var) const {return "update" + std::to_string(i_var);}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

}
#endif
