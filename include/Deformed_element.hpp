#ifndef DEFORMED_ELEMENT_HPP_
#define DEFORMED_ELEMENT_HPP_

#include <Element.hpp>

namespace cartdg
{

class Deformed_element : public Element
{
  int n_qpoint;
  Eigen::VectorXd jac;

  public:
  Deformed_element(Storage_params);
  double* jacobian();
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
};

class Deformed_elem_con
{
  public:
  std::array<Deformed_element*, 2> elems;
  std::array<int, 2> i_dims;
};

typedef std::vector<std::unique_ptr<Deformed_element>> def_elem_vec;
typedef std::vector<Deformed_elem_con> def_elem_con_vec;

}
#endif
