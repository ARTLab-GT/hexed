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

}
#endif
