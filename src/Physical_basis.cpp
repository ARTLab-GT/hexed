#include <limits>
#include <cmath>
#include <numeric>
#include <math.hpp>
#include <Physical_basis.hpp>

namespace cartdg
{

Physical_basis::Physical_basis(int n_dim, int row_size, std::vector<double> qpoint_pos)
: n_dim(n_dim), row_size(row_size), n_qpoint(custom_math::pow(row_size, n_dim)), pos(qpoint_pos),
  bound_box_center(n_dim), bound_box_size(n_dim), indices(n_dim)
{
  if (n_dim > 3) throw std::runtime_error("Demand for `Physical_basis` with `n_dim > 3` which is not supported.");
  if (int(pos.size()) != n_dim*n_qpoint) throw std::runtime_error("Size of `qpoint_pos` does not match `n_dim` and `row_size` in `Physical_basis.`");

  // compute bounding box
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    double min =  std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::max();
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      min = std::min(min, pos[i_dim*n_qpoint + i_qpoint]);
      max = std::max(max, pos[i_dim*n_qpoint + i_qpoint]);
    }
    bound_box_center[i_dim] = (min + max)/2.;
    bound_box_size[i_dim] = min + max;
  }

  // compute indices of polynomials
  std::vector<int> inds (n_dim, 0);
  for (int i_basis = 0; i_basis < size(); ++i_basis) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) indices[i_dim].push_back(inds[i_dim]);
    ++inds.back();
    for (int i_dim = n_dim - 1; i_dim > 0; --i_dim) {
      if (std::accumulate(inds.begin(), inds.end(), 0) >= row_size) {
        inds[i_dim] = 0;
        ++inds[i_dim - 1];
      }
    }
  }
}

int Physical_basis::size()
{
  switch (n_dim) {
    case 1: return row_size; break;
    case 2: return row_size*(row_size + 1)/2; break;
    case 3: return row_size*(row_size + 1)*(row_size + 2)/6; break;
    default: return 0;
  }
}

double Physical_basis::evaluate(int i_qpoint, int i_basis)
{
  double product = 1.;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    // scale/translate so that coordinates are +- 1 at edges of bounding box
    double transformed = (pos[i_dim*n_qpoint + i_qpoint] - bound_box_center[i_dim])*2/bound_box_size[i_dim];
    product *= std::legendre(indices[i_dim][i_basis], transformed);
  }
  return product;
}

Eigen::MatrixXd Physical_basis::projection(Eigen::MatrixXd polys, Eigen::VectorXd weights)
{
  return Eigen::MatrixXd::Zero(polys.rows(), polys.cols());
}

}
