#ifndef HEXED_ELEMENT_FUNC_HPP_
#define HEXED_ELEMENT_FUNC_HPP_

#include "Output_data.hpp"
#include "Qpoint_func.hpp"

namespace hexed
{

class Element_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const = 0;
};

// compute the average of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_average : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_average(Qpoint_func& func);
  virtual int n_var(int n_dim) const;
  virtual inline std::string variable_name(int i_var) const {return qf.variable_name(i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

// compute the L2 norm of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_l2 : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_l2(Qpoint_func&);
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual inline std::string variable_name(int i_var) const {return qf.variable_name(i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

}
#endif
