#include "Grid.hpp"

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_elem(n_elem_arg), basis(basis_arg)
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.rank;
  }
  n_dof = n_qpoint*n_var;
  storage_r.resize(n_dof*n_elem, 0.);
  access_r = storage_r.data();
  storage_w.resize(n_dof*n_elem, 0.);
  access_w = storage_w.data();
}

Grid::~Grid() {}
