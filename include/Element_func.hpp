#ifndef HEXED_ELEMENT_FUNC_HPP_
#define HEXED_ELEMENT_FUNC_HPP_

#include "Output_data.hpp"
#include "Qpoint_func.hpp"

namespace hexed
{

//! a function that has a single value for each element
class Element_func : virtual public Qpoint_func
{
  std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const override;
  public:
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const = 0; //!< \details --
};

/*! \brief function to fetch the value of the `resolution_badness` member of the `Element`.
 * \details Only useful if you have already set the `resolution_badness` member to something you are interested in.
 */
class Resolution_badness : virtual public Element_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const
  {
    return {elem.resolution_badness};
  }
};

//! computes the average of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_average : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_average(const Qpoint_func& func); //!< \param func function you want to compute the average of
  Elem_average(Qpoint_func&&) = delete; //!< can't accept temporaries because that could create a dangling reference
  virtual int n_var(int n_dim) const;
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "average_" + qf.variable_name(n_dim, i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

//! computes the \f$L_2\f$ norm of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_l2 : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_l2(const Qpoint_func&); //!< \param func function you want to compute the norm of
  Elem_l2(Qpoint_func&&) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "l2_" + qf.variable_name(n_dim, i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

/*! \brief compute the elementwise nonsmoothness indicator of the provided function
 * \details The elementwise nonsmoothness indicator
 * (different from the nonsmoothness indicator in `Solver::set_art_visc_smoothness`)
 * is essentially a measure of how much of the flow state
 * resides in the highest-order polynomial components.
 * If you want the nonsmoothness of the variable \f$u\f$, it is given by
 * \f[
 *    \sqrt{\frac{1}{n_{dim}} \sum_{0 \le i_{dim} < n_{dim}} \| \langle u(\vec{x}), P_{r - 1}(x_{i_{dim}}) \rangle P_{r - 1}(x_{i_{dim}}) \|^2}
 * \f]
 * where \f$P_j\f$ is the \f$j\f$th Legendre polynomial and \f$r\f$ is the row size of the basis.
 * In other words, you project the solution onto the highest-order univariate Legendre polynomial along each dimension,
 * take the \f$L_2\f$ norm of this in the other dimensions, and then take the RMS of that whole expression over all dimensions.
 */
class Elem_nonsmooth : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_nonsmooth(const Qpoint_func& func); //!< \param func the function you want to compute the nonsmoothness of (\f$u\f$ in the explanation above).
  Elem_nonsmooth(Qpoint_func&&) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual inline std::string variable_name(int n_dim, int i_var) const {return "nonsmoothness_" + qf.variable_name(n_dim, i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;

};

}
#endif
