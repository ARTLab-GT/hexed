#ifndef CARTDG_NEIGHBOR_AV_HPP_
#define CARTDG_NEIGHBOR_AV_HPP_

#include "jump.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void neighbor_av(double** connections_r, double** connections_w, int n_con,
                 int i_var, int i_axis,
                 const Eigen::VectorXd weights_1d, Kernel_settings& settings)
{
  jump<1, n_var, n_qpoint, row_size>(connections_r, connections_w, n_con,
                                     0, i_var, i_axis, weights_1d, settings);
}

}
#endif
