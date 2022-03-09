#include <iostream>

#include <Grid.hpp>
#include <math.hpp>
#include <Tecplot_file.hpp>

namespace cartdg
{

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_vertices(custom_math::pow(2, n_dim)), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), storage_params{3, n_var, n_dim, basis.row_size}
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    n_qpoint *= basis.row_size;
  }
  n_dof = n_qpoint*n_var;
  pos.resize(n_elem*n_dim, 0);
}

// Note: the return value is trivial right now, since it is equal to n_elem.
// However, once we implement adaptive refinement, elements may be created at places
// other than the end of the vector, and the return value will become non-trivial.
// Treat the return value as a black box for forward compatability.
int Grid::add_element(std::vector<int> position)
{
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) pos.push_back(position[i_dim]);
  return n_elem++;
}

void Grid::match_hanging(Element::shareable_value_access access_func)
{
  for (auto& matcher : hanging_matchers) {
    matcher.match(access_func);
  }
}

void Grid::print()
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> pos = get_pos(i_elem);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        std::cout << pos[i_qpoint + n_qpoint*i_dim] << "  ";
      }
    }
    std::cout << '\n';
  }
}

}
