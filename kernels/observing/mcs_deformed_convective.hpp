#ifndef CARTDG_MCS_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_MCS_DEFORMED_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Deformed_element.hpp>
#include "char_speed_convective.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 1)
template<int n_var, int n_qpoint, int row_size>
double mcs_deformed_convective(def_elem_vec& def_elements, Kernel_settings& settings)
{
  const int n_dim = n_var - 2;
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  double max_speed = 0.;
  #pragma omp parallel for reduction(max:max_speed)
  for (unsigned i_elem = 0; i_elem < def_elements.size(); ++i_elem)
  {
    double speed = char_speed_convective<n_var, n_qpoint>(def_elements[i_elem]->stage(i_read), heat_rat);
    double* jac_data = def_elements[i_elem]->jacobian();
    // account for jacobian
    double min_sv = 1.;
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      Eigen::Matrix<double, n_dim, n_dim> jac_mat;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          jac_mat(i_dim, j_dim) = jac_data[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
        }
      }
      // find minimum singular value (always at end of singular value vector)
      min_sv = std::min(min_sv, jac_mat.jacobiSvd().singularValues()(n_dim - 1));
    }
    // for "max speed" use maximum actual characteristic speed divided by minimum grid size,
    // the latter measured by the minimum singular value of the jacobian in the element
    max_speed = std::max<double>(speed/min_sv, max_speed);
  }
  return max_speed;
}

}
#endif
