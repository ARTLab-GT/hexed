#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include "Grid.hpp"
#include "Element.hpp"

namespace cartdg
{

class Qpoint_func
{
  public:
  virtual std::vector<double> operator()(Element& element, int i_qpoint) = 0;
};

class Jacobian_det_func
{
  public:
  virtual std::vector<double> operator()(Element& element, int i_qpoint);
};

}
#endif
