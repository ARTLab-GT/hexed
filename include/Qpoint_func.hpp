#ifndef HEXED_QPOINT_FUNC_HPP_
#define HEXED_QPOINT_FUNC_HPP_

#include <vector>
#include "Output_data.hpp"
#include "Element.hpp"
#include "Basis.hpp"

namespace hexed
{

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

// Determinant of inverse of Jacobian
class Jac_inv_det_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "jacobian_inverse_determinant";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

class Time_step_scale_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "time_step_scale";}
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

class Art_visc_coef : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "artificial_viscosity_coefficient";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

class Component : public Qpoint_func
{
  const Qpoint_func& qf;
  int iv;

  public:
  inline Component(const Qpoint_func& base, int i_var) : qf{base}, iv{i_var} {}
  inline Component(Qpoint_func&& base, int i_var) = delete; // can't accept temporaries. would cause dangling reference
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return qf.variable_name(iv);}
  virtual inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const
  {
    return {qf(e, b, i_qpoint, time)[iv]};
  }
};

class Scaled : public Qpoint_func
{
  const Qpoint_func& qf;
  std::array<double, 2> bnd;

  public:
  inline Scaled(const Qpoint_func& base, std::array<double, 2> bounds) : qf{base}, bnd{bounds} {}
  Scaled(Qpoint_func&& base, std::array<double, 2> bounds) = delete;
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int i_var) const {return "scaled_" + qf.variable_name(0);}
  virtual inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const
  {
    double val = qf(e, b, i_qpoint, time)[0];
    val = (val - bnd[0])/(bnd[1] - bnd[0]);
    return {val};
  }
};

class Pow : public Qpoint_func
{
  const Qpoint_func& qf;
  int exp;

  public:
  inline Pow(const Qpoint_func& base, int exponent) : qf{base}, exp{exponent} {}
  Pow(Qpoint_func&&, int) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual std::string variable_name(int i_var) const;
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

}
#endif
