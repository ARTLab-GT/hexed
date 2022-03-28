#include <Solver.hpp>

namespace cartdg
{

Solver::Solver(int n_dim, int row_size, double root_mesh_size) :
  params{3, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size}
{}

void Solver::initialize(const Spacetime_func& func)
{
}

std::vector<double> Solver::integral_field(const Qpoint_func& integrand)
{
  return {};
}

}
