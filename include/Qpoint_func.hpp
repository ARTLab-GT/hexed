#ifndef CARTDG_QPOINT_FUNC_HPP_
#define CARTDG_QPOINT_FUNC_HPP_

#include <vector>
#include "Output_data.hpp"

namespace cartdg
{

class Grid;

/*
 * Represents a function which can be evaluated at quadrature points. Can depend on
 * flow state, position, time, or on mathematical parameters like element Jabobian,
 * quadrature weights, etc. Note: `Output_data` is a virtual base class because
 * derived classes may inherit from multiple types of `Output_data`.
 */
class Qpoint_func : virtual public Output_data
{
  public:
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint) = 0;
};

// Returns a vector with one element: the Jacobian determinant at the quadrature point.
class Jacobian_det_func : public Qpoint_func
{
  public:
  virtual constexpr int n_var(int n_dim) {return 1;}
  virtual inline std::string variable_name(int i_var) {return "jacobian_determinant";}
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint);
};

class Time_step_scale_func : public Qpoint_func
{
  public:
  virtual constexpr int n_var(int n_dim) {return 1;}
  virtual inline std::string variable_name(int i_var) {return "time_step_scale";}
  virtual std::vector<double> operator()(Grid& grid, int i_element, int i_qpoint);
};

}
#endif
