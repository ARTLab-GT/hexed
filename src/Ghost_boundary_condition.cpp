#include <Ghost_boundary_condition.hpp>
#include <Grid.hpp>
#include <iostream>

namespace cartdg
{

Ghost_boundary_condition::Ghost_boundary_condition(const Grid& grid, int i_dim_arg,
                                                     bool is_positive_face_arg)
: i_dim(i_dim_arg), n_var(grid.n_var), n_qpoint(grid.n_qpoint/grid.basis.rank),
  is_positive_face(is_positive_face_arg), state(n_qpoint, 2*n_var)
{
  default_jacobian.clear();
  default_jacobian.resize(grid.n_dim*grid.n_dim*grid.n_qpoint, 0.);
  for (int i_dim = 0; i_dim < grid.n_dim; ++i_dim)
  {
    for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
    {
      default_jacobian[i_dim*(grid.n_dim + 1)*grid.n_qpoint + i_qpoint] = 1.;
    }
  }
}

void Ghost_boundary_condition::add_element(int i_elem, double* jacobian)
{
  elems.push_back(i_elem);
  jacobians.push_back(jacobian);
}

void Ghost_boundary_condition::add_element(int i_elem)
{
  add_element(i_elem, default_jacobian.data());
}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::domain_state()
{
  return (is_positive_face) ? state.block(0,     0, n_qpoint, n_var)
                            : state.block(0, n_var, n_qpoint, n_var);
}

Eigen::Block<Eigen::ArrayXXd> Ghost_boundary_condition::ghost_state()
{
  return (is_positive_face) ? state.block(0, n_var, n_qpoint, n_var)
                            : state.block(0,     0, n_qpoint, n_var);
}

void Ghost_boundary_condition::print()
{
  std::cout << i_dim << " " << is_positive_face << " : ";
  for (int i_elem : elems) std::cout << i_elem << " ";
  std::cout << "\n";
}

}
