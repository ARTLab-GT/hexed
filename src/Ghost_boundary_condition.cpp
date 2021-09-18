#include <Ghost_boundary_condition.hpp>
#include <Grid.hpp>
#include <iostream>

namespace cartdg
{

Ghost_boundary_condition::Ghost_boundary_condition(const Grid& grid, int i_dim_arg,
                                                     bool is_positive_face_arg)
: i_dim(i_dim_arg), n_var(grid.n_var), n_qpoint(grid.n_qpoint/grid.basis.row_size),
  is_positive_face(is_positive_face_arg), state(n_qpoint, 2*n_var)
{}

void Ghost_boundary_condition::add_element(int i_elem)
{
  elems.push_back(i_elem);
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
