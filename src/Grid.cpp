#include <iostream>

#include <Grid.hpp>
#include <math.hpp>
#include <get_mcs_cpg_euler.hpp>

namespace cartdg
{

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_vertices(custom_math::pow(2, n_dim)), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), i_rk_stage(0), i_read(0), i_write(1),
storage_params{3, n_var, n_dim, basis.row_size}
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.row_size;
    origin.push_back(0.);
  }
  n_dof = n_qpoint*n_var;
  for (int i = 0; i < (int)state_storage.size(); ++i) state_storage[i].resize(n_dof*n_elem, 0.);
  derivs.resize(n_elem*n_qpoint, 0.);
  visc.resize(n_elem*n_vertices, 0.);
  pos.resize(n_elem*n_dim, 0);
  for (int i_dim = 0; i_dim < 3*n_dim; ++i_dim)
  {
    neighbor_storage.emplace_back();
    deriv_neighbor_storage.emplace_back();
    visc_neighbor_storage.emplace_back();
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    viscous_inds.push_back(i_elem);
  }
}

double* Grid::state_r()
{
  return state_storage[i_read].data();
}

double* Grid::state_w()
{
  return state_storage[i_write].data();
}

void Grid::clear_neighbors()
{
  for (int i_dim = 0; i_dim < 2*n_dim; ++i_dim)
  {
    neighbor_storage[i_dim].clear();
    deriv_neighbor_storage[i_dim].clear();
    visc_neighbor_storage[i_dim].clear();
  }
}

// Note: the return value is trivial right now, since it is equal to n_elem.
// However, once we implement adaptive refinement, elements may be created at places
// other than the end of the vector, and the return value will become non-trivial.
// Treat the return value as a black box for forward compatability.
int Grid::add_element(std::vector<int> position)
{
  for (int i_rk_stage = 0; i_rk_stage < 3; ++i_rk_stage)
  {
    state_storage[i_rk_stage].resize(state_storage[i_rk_stage].size() + n_dof, 0.);
  }
  derivs.resize(derivs.size() + n_qpoint, 0.);
  visc.resize(visc.size() + n_vertices, 0.);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) pos.push_back(position[i_dim]);
  return n_elem++;
}

double Grid::stable_time_step(double cfl_by_stable_cfl, Kernel_settings& settings)
{
  double cfl = cfl_by_stable_cfl*get_stable_cfl();
  settings.i_read = i_read;
  return cfl*mesh_size/get_mcs_cpg_euler(n_dim, basis.row_size)(state_r(), n_elem, settings);
}

bool Grid::execute_runge_kutta_stage()
{
  double* read0 = state_storage[0].data();
  double* read1 = state_w();
  double* write = (i_rk_stage == 2) ? read0 : read1;
  const int r0 = 0;
  const int r1 = i_write;
  const int w = (i_rk_stage == 2) ? r0 : r1;
  double weight1 = rk_weights[i_rk_stage]; double weight0 = 1. - weight1;
  #pragma omp parallel for
  for (int i = 0; i < n_elem*n_dof; ++i)
  {
    write[i] = weight1*read1[i] + weight0*read0[i];
  }
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
  if ((basis.row_size > 0) && (basis.row_size <= 9))
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
      std::cout << ":    ";
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        std::cout << state_r()[i_qpoint + n_qpoint*i_var + n_dof*i_elem] << "    ";
      }
      std::cout << '\n';
    }
    std::cout << '\n';
  }
}

std::vector<double> Grid::integral()
{
  State_variables state_variables;
  return integral(state_variables);
}
std::vector<double> Grid::integral(Domain_func& integrand)
{
  Eigen::VectorXd weights (n_qpoint);
  Eigen::VectorXd weights_1d = basis.node_weights();
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) weights(i_qpoint) = 1.;
  for (int stride = n_qpoint/basis.row_size, n_rows = 1; n_rows < n_qpoint;
       stride /= basis.row_size, n_rows *= basis.row_size)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint)
        {
          weights((i_outer*basis.row_size + i_qpoint)*stride + i_inner)
          *= weights_1d(i_qpoint)*mesh_size;
        }
      }
    }
  }

  std::vector<double> total;
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> elem_pos = get_pos(i_elem);
    double* stage = element(i_elem).stage(0);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      std::vector<double> point_pos;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        point_pos.push_back(elem_pos[i_qpoint + i_dim*n_qpoint]);
      }
      std::vector<double> point_state;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        point_state.push_back(stage[i_qpoint + i_var*n_qpoint]);
      }
      std::vector<double> point_integrand = integrand(point_pos, time, point_state);
      double jac_det = jacobian_det(i_elem, i_qpoint);
      if (total.size() < point_integrand.size())
      {
        total.resize(point_integrand.size());
      }
      for (int i_var = 0; i_var < int(point_integrand.size()); ++i_var)
      {
        total[i_var] += point_integrand[i_var]*weights[i_qpoint]*jac_det;
      }
    }
  }
  return total;
}

}
