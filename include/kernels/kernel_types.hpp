#ifndef KERNEL_TYPES_HPP_
#define KERNEL_TYPES_HPP_

#include <Eigen/Dense>

typedef void (*Local_kernel)(double*, double*, int, const Eigen::MatrixXd&, double, double);
typedef void (*Neighbor_kernel)(double***, double***, int*, const Eigen::VectorXd,
              double, double);

#endif
