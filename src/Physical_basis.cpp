#include <math.hpp>
#include <Physical_basis.hpp>

namespace cartdg
{

Physical_basis::Physical_basis(int n_dim, int row_size, std::vector<double> qpoint_pos)
: n_dim(n_dim), row_size(row_size), n_qpoint(custom_math::pow(row_size, n_dim)), pos(qpoint_pos)
{
  if (n_dim > 3) throw std::runtime_error("Demand for `Physical_basis` with `n_dim > 3` which is not supported.");
  if (int(pos.size()) != n_dim*n_qpoint) throw std::runtime_error("Size of `qpoint_pos` does not match `n_dim` and `row_size` in `Physical_basis.`");
}

}
