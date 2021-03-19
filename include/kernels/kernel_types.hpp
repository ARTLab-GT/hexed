#ifndef CARTDG_KERNEL_TYPES_HPP_
#define CARTDG_KERNEL_TYPES_HPP_

#include <Eigen/Dense>

#include "Kernel_settings.hpp"

namespace cartdg
{

typedef void (*Local_kernel)(double*, double*, int, Eigen::MatrixXd, Kernel_settings&);
typedef void (*Neighbor_kernel)(double***, double***, int*, const Eigen::VectorXd,
              Kernel_settings&);
typedef void (*Gbc_kernel)(std::vector<Ghost_boundary_condition*>&, double*, double*,
                           double, Kernel_settings&);
typedef void (*Nonpen_kernel)(double*, double*, double*, int*, int*, int*, int, double,
                              Kernel_settings&);
typedef double (*Max_char_speed_kernel)(double*, int, Kernel_settings&);
typedef double (*Physical_step_kernel)(double*, double*, int, Kernel_settings&);
typedef void (*Restrict_step_kernel)(double*, double*, int, double, Kernel_settings&);
typedef void (*Local_deformed_kernel)(double*, double*, double*, int,
                                      const Eigen::MatrixXd&, Kernel_settings&);
typedef void (*Neighbor_deformed_kernel)(double**, double**, double**, int*, int*,
                                         int, const Eigen::VectorXd, Kernel_settings&);
                             

}
#endif
