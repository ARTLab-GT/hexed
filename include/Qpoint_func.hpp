#ifndef HEXED_QPOINT_FUNC_HPP_
#define HEXED_QPOINT_FUNC_HPP_

#include <vector>
#include "Output_data.hpp"
#include "Element.hpp"
#include "Basis.hpp"

namespace hexed
{

/*! \brief Represents a function which can be evaluated at quadrature points.
 * \details Can depend on flow state, position, time, or on mathematical parameters like element Jabobian,
 * quadrature weights, etc. Note: `Output_data` is a virtual base class because
 * derived classes may inherit from multiple types of `Output_data`.
 */
class Qpoint_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const = 0;
};

//! Returns a vector with one element: the Jacobian determinant at the quadrature point.
class Jacobian_det_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "jacobian_determinant";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! Determinant of inverse of Jacobian
class Jac_inv_det_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "jacobian_inverse_determinant";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! Fetches the value of the `Element::time_step_scale` member.
class Time_step_scale_func : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "time_step_scale";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! returns the most recent update to the state divided by the local time step scale.
class Physical_update : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return n_dim + 2;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "update" + std::to_string(i_var);}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! fetches the artificial viscosity coefficient
class Art_visc_coef : public Qpoint_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "artificial_viscosity_coefficient";}
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! returns a component of another `Qpoint_func`
class Component : public Qpoint_func
{
  const Qpoint_func& qf;
  int iv;

  public:
  /*! \param base Function you want to get a component of
   * \param i_var Which component you want
   */
  inline Component(const Qpoint_func& base, int i_var) : qf{base}, iv{i_var} {}
  inline Component(Qpoint_func&& base, int i_var) = delete; // can't accept temporaries. would cause dangling reference
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return qf.variable_name(n_dim, iv);}
  virtual inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const
  {
    return {qf(e, b, i_qpoint, time)[iv]};
  }
};

/*! \brief Scales a function of a single variable.
 * \details Computes a linear transformation that maps `bounds` to [0, 1]
 * and then applies this transformation to the input function (`base`).
 * Thus if `bounds` are lower and upper bounds for the function,
 * then `Scaled` will return an output which is in the range [0, 1].
 * This is useful for manually computing contour plots or colormapping data.
 */
class Scaled : public Qpoint_func
{
  const Qpoint_func& qf;
  std::array<double, 2> bnd;

  public:
  /*!
   * \param base Function you want to scale.
   * \param bounds Points to define the linear map (see class description).
   */
  inline Scaled(const Qpoint_func& base, std::array<double, 2> bounds) : qf{base}, bnd{bounds} {}
  Scaled(Qpoint_func&& base, std::array<double, 2> bounds) = delete;
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "scaled_" + qf.variable_name(n_dim, 0);}
  virtual inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const
  {
    double val = qf(e, b, i_qpoint, time)[0];
    val = (val - bnd[0])/(bnd[1] - bnd[0]);
    return {val};
  }
};

//! Raises the output of a `Qpoint_func` to a user-specified power.
class Pow : public Qpoint_func
{
  const Qpoint_func& qf;
  int exp;

  public:
  inline Pow(const Qpoint_func& base, int exponent) : qf{base}, exp{exponent} {}
  Pow(Qpoint_func&&, int) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual std::string variable_name(int n_dim, int i_var) const;
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
};

//! Concatenates `Qpoint_func`s
typedef Concat_func<Qpoint_func, Element&, const Basis&, int, double> Qf_concat;

}
#endif
