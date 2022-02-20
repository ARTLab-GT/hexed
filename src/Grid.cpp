#include <iostream>

#include <Grid.hpp>
#include <math.hpp>
#include <Tecplot_file.hpp>

namespace cartdg
{

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_vertices(custom_math::pow(2, n_dim)), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), i_rk_stage(0), i_read(0), i_write(1),
storage_params{3, n_var, n_dim, basis.row_size}
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    n_qpoint *= basis.row_size;
  }
  n_dof = n_qpoint*n_var;
  pos.resize(n_elem*n_dim, 0);
}

int Grid::i_stage_read()
{
  return i_read;
}

int Grid::i_stage_write()
{
  return i_write;
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

bool Grid::execute_runge_kutta_stage()
{
  const int r0 = 0;
  const int r1 = i_write;
  const int w = (i_rk_stage == 2) ? r0 : r1;
  double weight1 = rk_weights[i_rk_stage]; double weight0 = 1. - weight1;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    auto& elem = element(i_elem);
    double* stage_r0 = elem.stage(r0);
    double* stage_r1 = elem.stage(r1);
    double* stage_w  = elem.stage(w );
    for (int i_dof = 0; i_dof < storage_params.n_dof(); ++i_dof)
    {
      stage_w[i_dof] = weight1*stage_r1[i_dof] + weight0*stage_r0[i_dof];
    }
  }
  ++i_rk_stage; ++i_write;
  if (i_write == 3) i_write = 1;
  if (i_rk_stage == 3) { i_rk_stage = 0; i_write = 1; ++iter; }
  i_read = i_rk_stage;
  return i_rk_stage == 0;
}

double Grid::get_stable_cfl()
{
  if ((basis.row_size > 0) && (basis.row_size <= 10))
  {
    return stable_cfl[basis.row_size - 1];
  }
  else
  {
    throw std::runtime_error("Stable CFL number unknown for basis of desired row_size.");
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
