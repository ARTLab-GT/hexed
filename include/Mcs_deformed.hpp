#ifndef HEXED_MCS_DEFORMED_HPP_
#define HEXED_MCS_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "kernel_factory.hpp"
#include "Deformed_element.hpp"
#include "math.hpp"
#include "characteristic_speed.hpp"

namespace hexed
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
      double* ref_nrml = elem.reference_level_normals();
      double* jac_det = elem.jacobian_determinant();
      double min_scale = std::numeric_limits<double>::max();
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        double norm_sum = 0.;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          double norm_sq = 0.;
          for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
            double coef = ref_nrml[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
            norm_sq += coef*coef;
          }
          norm_sum += std::sqrt(norm_sq);
        }
        double determinant = jac_det[i_qpoint];
        min_scale = std::min(min_scale, n_dim*determinant/tss[i_qpoint]/norm_sum);
      }
      // compute reference speed
      double speed = characteristic_speed<n_dim, n_qpoint>(elem.stage(0), heat_rat)/min_scale/elem.nominal_size();
      max_speed = std::max(speed, max_speed);
    }
    return max_speed;
  }
};

}
#endif
