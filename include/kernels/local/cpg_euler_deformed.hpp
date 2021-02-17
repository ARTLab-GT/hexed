#ifndef CARTDG_CPG_EULER_DEFORMED_HPP_
#define CARTDG_CPG_EULER_DEFORMED_HPP_

#include <Eigen/Dense>
#include "../Kernel_settings.hpp"

namespace cartdg
{

template<int n_var, int n_qpoint, int row_size>
void cpg_euler_deformed(double* read, double* write, double* jacobian, int n_elem,
                        const Eigen::MatrixXd& diff_mat_arg,
                        Kernel_settings& settings)
{}

}

#endif