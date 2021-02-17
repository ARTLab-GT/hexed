#ifndef CARTDG_KERNEL_TYPES_HPP_
#define CARTDG_KERNEL_TYPES_HPP_

#include <Eigen/Dense>

#include "Kernel_settings.hpp"

namespace cartdg
{

typedef void (*Local_kernel)(double*, double*, int, const Eigen::MatrixXd&, Kernel_settings&);
typedef void (*Neighbor_kernel)(double***, double***, int*, const Eigen::VectorXd,
              Kernel_settings&);
typedef void (*Fbc_kernel)(std::vector<Fitted_boundary_condition*>&, double*, double*,
                           double, Kernel_settings&);
typedef double (*Max_char_speed_kernel)(double*, int, Kernel_settings&);

}
#endif
