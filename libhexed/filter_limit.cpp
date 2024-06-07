#include <filter_limit.hpp>
#include <math.hpp>
#include <Row_index.hpp>

namespace hexed
{

void filter_limit(int n_dim, double* data, const Basis& basis, double decay_rate)
{
  int row_size = basis.row_size;
  int n_qpoint = math::pow(row_size, n_dim);
  Eigen::Map<Mat<>> qpoints(data, n_qpoint);
  Mat<dyn, dyn> orthogonal(row_size, row_size);
  for (int i_row = 0; i_row < row_size; ++i_row) {
    orthogonal(i_row, all) = basis.orthogonal(i_row);
  }
  Mat<dyn, dyn> weighted = orthogonal*basis.node_weights().asDiagonal();
  Mat<> modes = math::hypercube_matvec(weighted, qpoints);
  double norm = modes.norm();
  for (int i_mode = 0; i_mode < n_qpoint; ++i_mode) {
    int total = 0;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) total += (i_mode/Row_index(n_dim, row_size, i_dim).stride)%row_size;
    double threshold = norm*std::pow(decay_rate, total);
    double mode_norm = std::abs(modes(i_mode)); // note each of the Legendre modes have mode_norm 1
    if (mode_norm > threshold) modes(i_mode) *= threshold/mode_norm;
  }
  qpoints = math::hypercube_matvec(orthogonal.transpose(), modes);
}

}
