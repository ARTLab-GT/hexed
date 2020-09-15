#include "Grid.hpp"

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg)
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.rank;
  }
  n_dof = n_qpoint*n_var;
  state_r.resize(n_dof*n_elem, 0.);
  state_w.resize(n_dof*n_elem, 0.);
  pos.resize(n_elem*n_dim, 0);
}

Grid::~Grid() {}

std::vector<double> Grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim, 0.);
  return elem_pos;
}

void Grid::visualize()
{
}
