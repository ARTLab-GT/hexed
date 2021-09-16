#include <Regular_grid.hpp>
#include <get_mcs_convective.hpp>
#include <get_local_convective.hpp>
#include <get_neighbor_convective.hpp>
#include <get_gbc_convective.hpp>
#include <get_req_visc_cpg_euler.hpp>
#include <get_cont_visc_cpg_euler.hpp>
#include <get_local_derivative.hpp>
#include <get_neighbor_derivative.hpp>
#include <get_av_flux.hpp>
#include <get_local_av.hpp>
#include <get_neighbor_av.hpp>
#include <get_gbc_av.hpp>

namespace cartdg
{

Regular_grid::Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: Grid (n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg), elements{}
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    elements.emplace_back(new Element {storage_params});
  }
  for (int i_dim = 0; i_dim < 3*n_dim; ++i_dim)
  {
    elem_cons.push_back({});
  }
}

Element& Regular_grid::element(int i_elem)
{
  return *elements[i_elem];
}

elem_con Regular_grid::connection(int i_dim, int i_con)
{
  return elem_cons[i_dim][i_con];
}

int Regular_grid::n_con(int i_dim)
{
  return int(elem_cons[i_dim].size());
}

std::vector<double> Regular_grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim, 0.);
  std::vector<int> indices;
  populate_slice(elem_pos, indices, i_elem);
  return elem_pos;
}

double Regular_grid::jacobian_det(int i_elem, int i_qpoint)
{
  return 1.;
}

std::vector<double**> Regular_grid::neighbor_connections_r()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(neighbor_storage[i_dim + i_read*n_dim].data());
  }
  return connections;
}

std::vector<double**> Regular_grid::deriv_neighbor_connections()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(deriv_neighbor_storage[i_dim].data());
  }
  return connections;
}

std::vector<double**> Regular_grid::visc_neighbor_connections()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(visc_neighbor_storage[i_dim].data());
  }
  return connections;
}

std::vector<double**> Regular_grid::neighbor_connections_w()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(neighbor_storage[i_dim + i_write*n_dim].data());
  }
  return connections;
}

std::vector<int> Regular_grid::n_neighb_con()
{
  std::vector<int> n_con;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_con.push_back(neighbor_storage[i_dim].size()/2);
  }
  return n_con;
}

double Regular_grid::stable_time_step(double cfl_by_stable_cfl, Kernel_settings& settings)
{
  double cfl = cfl_by_stable_cfl*get_stable_cfl();
  settings.i_read = i_read;
  return cfl*mesh_size/get_mcs_convective(n_dim, basis.row_size)(elements, settings);
}

void Regular_grid::execute_local(Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_local_convective(n_dim, basis.row_size)(elements, basis, settings);
}

void Regular_grid::execute_neighbor(Kernel_settings& settings)
{
  settings.i_read = i_read;
  settings.i_write = i_write;
  get_neighbor_convective(n_dim, basis.row_size)(elem_cons, basis, settings);
  //get_gbc_cpg_euler(n_dim, basis.row_size)(ghost_bound_conds, state_r(), state_w(), basis, settings);
}

void Regular_grid::execute_req_visc(Kernel_settings& settings)
{
  get_req_visc_cpg_euler(n_dim, basis.row_size)(state_r(), visc.data(), n_elem, basis, settings);
}

void Regular_grid::execute_cont_visc(Kernel_settings& settings)
{
  get_cont_visc_cpg_euler(n_dim, basis.row_size)(visc_neighbor_connections().data(), n_neighb_con().data(), settings);
}

void Regular_grid::execute_local_derivative(int i_var, int i_dim, Kernel_settings& settings)
{
  get_local_derivative(n_dim, basis.row_size)(state_r(), derivs.data(), n_elem, i_var, i_dim, basis, settings);
}

void Regular_grid::execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings& settings)
{
  int n_con = n_neighb_con()[i_dim];
  get_neighbor_derivative(n_dim, basis.row_size)(neighbor_connections_r()[i_dim], deriv_neighbor_connections()[i_dim], n_con, i_var, i_dim, basis, settings);
}

void Regular_grid::execute_av_flux(Kernel_settings& settings)
{
  get_av_flux(n_dim, basis.row_size)(derivs.data(), visc.data(), n_elem, basis, settings);
}

void Regular_grid::execute_local_av(int i_var, int i_dim, Kernel_settings& settings)
{
  get_local_av(n_dim, basis.row_size)(derivs.data(), state_w(), n_elem, i_var, i_dim, basis, settings);
}

void Regular_grid::execute_neighbor_av(int i_var, int i_dim, Kernel_settings& settings)
{
  int n_con = n_neighb_con()[i_dim];
  get_neighbor_av(n_dim, basis.row_size)(deriv_neighbor_connections()[i_dim],
                                         neighbor_connections_w()[i_dim], n_con, i_var, i_dim, basis, settings);
  get_gbc_av(n_dim, basis.row_size)(ghost_bound_conds, derivs.data(), state_w(), i_var, i_dim, basis, settings);
}

int Regular_grid::add_element(std::vector<int> position)
{
  elements.emplace_back(new Element {storage_params});
  return Grid::add_element(position);
}

void Regular_grid::add_connection(int i_elem0, int i_elem1, int i_dim)
{
  for (int i_stage = 0; i_stage < 3; ++i_stage)
  {
    double* state_data = state_storage[i_stage].data();
    for (int i_elem : {i_elem0, i_elem1})
    {
      neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data + n_dof*i_elem);
    }
  }
  for (int i_elem : {i_elem0, i_elem1})
  {
    deriv_neighbor_storage[i_dim].push_back(derivs.data() + n_qpoint*i_elem);
    visc_neighbor_storage[i_dim].push_back(visc.data() + n_vertices*i_elem);
  }
  elem_cons[i_dim].push_back({elements[i_elem0].get(), elements[i_elem1].get()});
}

void Regular_grid::clear_neighbors()
{
  Grid::clear_neighbors();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    elem_cons[i_dim].clear();
  }
}

void Regular_grid::auto_connect(std::vector<int> periods)
{
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    for (int j_elem = i_elem + 1; j_elem < n_elem; ++j_elem)
    {
      int pos_diff [n_dim];
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        pos_diff[i_dim] = pos[n_dim*j_elem + i_dim] - pos[n_dim*i_elem + i_dim];
        // assumes period[i_dim] >= greatest distance between elements in this dimension
        if (periods[i_dim] > 0)
        {
          pos_diff[i_dim] = (pos_diff[i_dim] + periods[i_dim] + 1)%periods[i_dim] - 1;
        }
      }
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        bool is_same_row = true;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim)
        {
          if ((j_dim != i_dim) && (pos_diff[j_dim] != 0)) is_same_row = false;
        }
        if (is_same_row)
        {
          if      (pos_diff[i_dim] ==  1) add_connection(i_elem, j_elem, i_dim);
          else if (pos_diff[i_dim] == -1) add_connection(j_elem, i_elem, i_dim);
        }
      }
    }
  }
}

void Regular_grid::auto_connect()
{
  std::vector<int> periods;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) periods.push_back(0);
  auto_connect(periods);
}

void Regular_grid::populate_slice(std::vector<double>& elem_pos, std::vector<int> indices, int i_elem)
{
  if ((int)indices.size() < n_dim)
  {
    indices.push_back(0);
    for (int i = 0; i < basis.row_size; ++i)
    {
      indices.back() = i;
      populate_slice(elem_pos, indices, i_elem);
    }
  }
  else
  {
    int i_flat = 0;
    int stride = n_qpoint;
    for (auto i : indices)
    {
      stride /= basis.row_size;
      i_flat += i*stride;
    }
    for (int i_dim = 0; i_dim < n_dim; ++i_dim)
    {
      elem_pos[i_flat + n_qpoint*i_dim] = (basis.node(indices[i_dim])
                                           + pos[i_elem*n_dim + i_dim])*mesh_size
                                           + origin[i_dim];
    }
  }
}

}
