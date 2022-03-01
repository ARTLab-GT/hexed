#ifndef CARTDG_REFERENCE_DERIVATIVE_HPP_
#define CARTDG_REFERENCE_DERIVATIVE_HPP_

#include <Eigen/Dense>

namespace cartdg
{

/*
 * Compute the gradient of a variable in reference space. That is, the data pointed to by `read` is treated as
 * though it describes a polynomial in the unit hypercube, with layout `[i_qpoint]`. All components of the gradient
 * are written into `write` with layout `[i_dim][i_qpoint]`. The differentiation is performed by applying `diff_mat`
 * to each hypercolumn.
 */
template<int n_dim, int n_qpoint, int row_size>
void reference_derivative(double* read, double* write, const Eigen::Matrix<double, row_size, row_size>& diff_mat)
{
}

}
#endif
