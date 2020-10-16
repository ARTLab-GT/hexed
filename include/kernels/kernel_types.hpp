#ifndef CARTDG_KERNEL_TYPES_HPP_
#define CARTDG_KERNEL_TYPES_HPP_

#include <Eigen/Dense>

namespace cartdg
{

typedef void (*Local_kernel)(double*, double*, int, const Eigen::MatrixXd&, double, double);
typedef void (*Neighbor_kernel)(double***, double***, int*, const Eigen::VectorXd,
                                double, double);
typedef void (*Ghost_function)(double* read, double* write, int n_points);
typedef void (*Neighbor_boundary_kernel)(double***, double***, int*, const Eigen::VectorXd,
                                         Ghost_function, double, double);
typedef double (*Max_char_speed_kernel)(double*, int, double);

}
#endif
