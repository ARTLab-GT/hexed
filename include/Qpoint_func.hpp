#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include <vector>

namespace cartdg
{

class Grid;

/*
 * Represents a function which can be evaluated at quadrature points. Can depend on
 * flow state, position, time, or on mathematical parameters like element Jabobian,
 * quadrature weights, etc.
 */
class Qpoint_func
{
  public:
  virtual ~Qpoint_func() = default;
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint) = 0;
};

// Returns a vector with one element: the Jacobian determinant at the quadrature point.
class Jacobian_det_func : public Qpoint_func
{
  public:
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint);
};

}
#endif
