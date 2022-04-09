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
template <int n_dim, int row_size>
class Mcs_deformed : public Kernel<Deformed_element&, double>
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
      // account for jacobian and time step scale
      double* tss = elem.time_step_scale();
      double* jac_data = elem.jacobian();
      double min_scale = std::numeric_limits<double>::max();
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        Eigen::Matrix<double, n_dim, n_dim> jac_mat;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
            jac_mat(i_dim, j_dim) = jac_data[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
          }
        }
        // find minimum singular value (always at end of singular value vector) and divide by TSS
        min_scale = std::min(min_scale, jac_mat.jacobiSvd().singularValues()(n_dim - 1)/tss[i_qpoint]);
      }
      // compute reference speed
      double speed = characteristic_speed<n_dim, n_qpoint>(elem.stage(0), heat_rat)/min_scale/elem.nominal_size();
      max_speed = std::max(speed, max_speed);
    }
    return max_speed;
  }
};

template<>
class Kernel_traits<Mcs_deformed>
{
  public:
  using base_t = Kernel<Deformed_element&, double>;
};

}
#endif
