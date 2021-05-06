#ifndef CARTDG_NEIGHBOR_DERIVATIVE_HPP_
#define CARTDG_NEIGHBOR_DERIVATIVE_HPP_

#include "jump.hpp"

namespace cartdg
{

// AUTOGENERATE
template<int n_var, int n_qpoint, int row_size>
void neighbor_derivative(double** connections_r, double** connections_w, int n_con,
                         int i_var, int i_axis,
                         const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  jump<n_var, 1, n_qpoint, row_size>(connections_r, connections_w, n_con,
                                     i_var, 0, i_axis, weights_1d, settings);
}

}
#endif
