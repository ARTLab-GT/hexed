#ifndef CARTDG_MCS_DEFORMED_HPP_
#define CARTDG_MCS_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "kernel_factory.hpp"
#include "Deformed_element.hpp"
#include "math.hpp"
#include "characteristic_speed.hpp"

namespace cartdg
{

/*
 * Computes the maximum speed of characteristics in reference space, accounting for any local time step.
 * For example, if the flow velocity is 0, the speed of sound is 340, the mesh size is 0.1,
 * the minimum singular value of the Jacobian is 0.5, and the local time step scale is 0.3,
 * then the computed characteristic speed will be 340*0.3/0.1/0.5 = 2040.
 */
class Mcs_deformed_dynamic
{
  public:
  virtual ~Mcs_deformed_dynamic() = default;
  virtual double operator()(Sequence<Deformed_element&>&) = 0;
};

template <int n_dim, int row_size>
class Mcs_deformed : public Mcs_deformed_dynamic
{
  double heat_rat;

  public:
  Mcs_deformed(double heat_ratio = 1.4) : heat_rat{heat_ratio} {}

  virtual double operator()(Sequence<Deformed_element&>& elements)
  {
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    double max_speed = 0.;
    #pragma omp parallel for reduction(max:max_speed)
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      Deformed_element& elem {elements[i_elem]};
      double max_tss = 0.;
      double* tss = elem.time_step_scale();
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        max_tss = std::max(max_tss, tss[i_qpoint]);
      }
      double speed = characteristic_speed<n_dim, n_qpoint>(elem.stage(0), heat_rat)*max_tss/elem.nominal_size();
      max_speed = std::max(speed, max_speed);
    }
    return max_speed;
  }
};

template<>
class Kernel_traits<Mcs_deformed>
{
  public:
  using base_t = Mcs_deformed_dynamic;
};

}
#endif
