#ifndef HEXED_ELEMENT_FUNC_HPP_
#define HEXED_ELEMENT_FUNC_HPP_

#include "Output_data.hpp"
#include "Qpoint_func.hpp"

namespace hexed
{

class Element_func : virtual public Qpoint_func
{
  virtual std::vector<double> operator()(Element&, const Basis&, int i_qpoint, double time) const;
  public:
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const = 0;
};

class Resolution_badness : virtual public Element_func
{
  public:
  virtual inline int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const
  {
    return {elem.resolution_badness};
  }
};

// compute the average of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_average : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_average(const Qpoint_func&);
  Elem_average(Qpoint_func&&) = delete; // can't accept temporaries because that could create a dangling reference
  virtual int n_var(int n_dim) const;
  virtual inline std::string variable_name(int i_var) const {return "average_" + qf.variable_name(i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

// compute the L2 norm of the provided `Qpoint_func` within the element by Gaussian quadrature
class Elem_l2 : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_l2(const Qpoint_func&);
  Elem_l2(Qpoint_func&&) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual inline std::string variable_name(int i_var) const {return "l2_" + qf.variable_name(i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;
};

// compute the elementwise nonsmoothness indicator of the provided function
class Elem_nonsmooth : public Element_func
{
  const Qpoint_func& qf;
  public:
  Elem_nonsmooth(const Qpoint_func&);
  Elem_nonsmooth(Qpoint_func&&) = delete;
  virtual inline int n_var(int n_dim) const {return qf.n_var(n_dim);}
  virtual inline std::string variable_name(int i_var) const {return "nonsmoothness_" + qf.variable_name(i_var);}
  virtual std::vector<double> operator()(Element& elem, const Basis&, double time) const;

};

}
#endif
