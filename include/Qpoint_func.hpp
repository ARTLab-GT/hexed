#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include "Grid.hpp"

namespace cartdg
{

class Qpoint_func
{
  protected:
  Grid& grid;
  public:
  Qpoint_func(Grid&);
  virtual std::vector<double> operator()(int i_element, int i_qpoint) = 0;
};

class Jacobian_det_func : public Qpoint_func
{
  public:
  Jacobian_det_func(Grid&);
  virtual std::vector<double> operator()(int i_element, int i_qpoint);
};

}
#endif
