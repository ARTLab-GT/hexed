#ifndef CARTDG_KERNEL_TYPES_HPP_
#define CARTDG_KERNEL_TYPES_HPP_

#include <Eigen/Dense>

#include "../Basis.hpp"
#include "Kernel_settings.hpp"

namespace cartdg
{

typedef void (*Nonpen_kernel)(double*, double*, double*, int*, int*, int*, int, double,
                              Kernel_settings&);
typedef void (*Derivative_kernel)(double*, double*, int, int, int, Basis&, Kernel_settings&);
typedef void (*Viscous_local_kernel)(double*, double*, int, int, int, Basis&, Kernel_settings&);
typedef void (*Jump_kernel)(double**, double**, int, int, int, Eigen::VectorXd, Kernel_settings&);
typedef void (*Jump_gbc_kernel)(std::vector<Ghost_boundary_condition*>&, double*, double*, int, int, double, Kernel_settings&);
typedef void (*Viscous_neighbor_kernel)(double**, double**, int, int, int, Eigen::VectorXd, Kernel_settings&);

}
#endif
