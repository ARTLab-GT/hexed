#include <limits>
#include <cmath>
#include <numeric>
#include <math.hpp>
#include <Physical_basis.hpp>

namespace hexed
{

Physical_basis::Physical_basis(int n_dim, int row_size, std::vector<double> qpoint_pos)
: n_dim(n_dim), row_size(row_size), n_qpoint(math::pow(row_size, n_dim)), pos(qpoint_pos),
  bound_box_center(n_dim), bound_box_size(n_dim), indices(n_dim)
{
  if (n_dim > 3) throw std::runtime_error("Demand for `Physical_basis` with `n_dim > 3` which is not supported.");
  if (int(pos.size()) != n_dim*n_qpoint) throw std::runtime_error("Size of `qpoint_pos` does not match `n_dim` and `row_size` in `Physical_basis.`");

  // compute bounding box
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    double max = *std::max_element(pos.begin() + i_dim*n_qpoint, pos.begin() + (i_dim + 1)*n_qpoint);
    double min = *std::min_element(pos.begin() + i_dim*n_qpoint, pos.begin() + (i_dim + 1)*n_qpoint);
    bound_box_center[i_dim] = (min + max)/2.;
    bound_box_size[i_dim] = max - min;
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
    double transformed = (pos[i_dim*n_qpoint + i_qpoint] - bound_box_center[i_dim])/bound_box_size[i_dim]*1.99; // multiply by slightly less than 2. to guard against rounding error
    product *= std::legendre(indices[i_dim][i_basis], transformed);
  }
  return product;
}

Eigen::MatrixXd Physical_basis::projection(Eigen::MatrixXd polys, Eigen::VectorXd weights)
{
  Eigen::MatrixXd mat (n_qpoint, size());
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
    for (int i_basis = 0; i_basis < size(); ++i_basis) {
      mat(i_qpoint, i_basis) = evaluate(i_qpoint, i_basis);
    }
  }
  auto sqrt_w = weights.cwiseSqrt();
  return mat*(sqrt_w.asDiagonal()*mat).fullPivHouseholderQr().solve(sqrt_w.asDiagonal()*polys);
}

}
