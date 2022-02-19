#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include <vector>

namespace cartdg
{

class Grid;

class Qpoint_func
{
  public:
  virtual ~Qpoint_func() = default;
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint) = 0;
};

class Jacobian_det_func : public Qpoint_func
{
  public:
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint);
};

}
#endif
