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
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "jacobian_determinant";}
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! Determinant of inverse of Jacobian
class Jac_inv_det_func : public Qpoint_func
{
  public:
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "jacobian_inverse_determinant";}
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! Fetches the value of the `Element::time_step_scale` member.
class Time_step_scale_func : public Qpoint_func
{
  public:
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "time_step_scale";}
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! returns the most recent update to the state divided by the local time step scale.
class Physical_update : public Qpoint_func
{
  public:
  inline int n_var(int n_dim) const override {return n_dim + 2;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "update" + std::to_string(i_var);}
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! fetches the artificial viscosity coefficient
class Art_visc_coef : public Qpoint_func
{
  public:
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "artificial_viscosity_coefficient";}
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
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
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return qf.variable_name(n_dim, iv);}
  inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const override
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
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "scaled_" + qf.variable_name(n_dim, 0);}
  inline std::vector<double> operator()(Element& e, const Basis& b, int i_qpoint, double time) const override
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
  inline int n_var(int n_dim) const override {return qf.n_var(n_dim);}
  std::string variable_name(int n_dim, int i_var) const override;
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! Fetches first forcing variable in artificial viscosity computation
class Art_visc_forcing : public Qpoint_func
{
  int i_force;
  public:
  //! \param i_forcing which forcing term you want to see
  inline Art_visc_forcing(int i_forcing) : i_force{i_forcing} {}
  inline int n_var(int n_dim) const override {return 1;}
  inline std::string variable_name(int n_dim, int i_var) const override {return "art_visc_forcing";};
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
};

//! Concatenates `Qpoint_func`s
typedef Concat_func<Qpoint_func, Element&, const Basis&, int, double> Qf_concat;

}
#endif
