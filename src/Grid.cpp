#include <iostream>

#include <Grid.hpp>

namespace cartdg
{

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), i_rk_stage(0), i_read(0), i_write(1)
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.rank;
    origin.push_back(0.);
  }
  n_dof = n_qpoint*n_var;
  for (int i = 0; i < (int)state_storage.size(); ++i) state_storage[i].resize(n_dof*n_elem, 0.);
  derivs.resize(n_elem*n_qpoint, 0.);
  pos.resize(n_elem*n_dim, 0);
  for (int i_dim = 0; i_dim < 3*n_dim; ++i_dim)
  {
    neighbor_storage.emplace_back();
    deriv_neighbor_storage.emplace_back();
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    viscous_inds.push_back(i_elem);
  }
}

Grid::~Grid() {}

double* Grid::state_r()
{
  return state_storage[i_read].data();
}

double* Grid::state_w()
{
  return state_storage[i_write].data();
}

std::vector<double**> Grid::neighbor_connections_r()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(neighbor_storage[i_dim + i_read*n_dim].data());
  }
  return connections;
}

std::vector<double**> Grid::deriv_neighbor_connections()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(deriv_neighbor_storage[i_dim].data());
  }
  return connections;
}

std::vector<double**> Grid::neighbor_connections_w()
{
  std::vector<double**> connections;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    connections.push_back(neighbor_storage[i_dim + i_write*n_dim].data());
  }
  return connections;
}

std::vector<int> Grid::n_neighb_con()
{
  std::vector<int> n_con;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_con.push_back(neighbor_storage[i_dim].size()/2);
  }
  return n_con;
}

void Grid::auto_connect(std::vector<int> periods)
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

void Grid::auto_connect()
{
  std::vector<int> periods;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) periods.push_back(0);
  auto_connect(periods);
}

void Grid::clear_neighbors()
{
  for (int i_dim = 0; i_dim < 2*n_dim; ++i_dim)
  {
    neighbor_storage[i_dim].clear();
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
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) pos.push_back(position[i_dim]);
  return n_elem++;
}

void Grid::add_connection(int i_elem0, int i_elem1, int i_dim)
{
  for (int i_stage = 0; i_stage < 3; ++i_stage)
  {
    double* state_data = state_storage[i_stage].data();
    for (int i_elem : {i_elem0, i_elem1})
    {
      neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data + n_dof*i_elem);
      deriv_neighbor_storage[i_dim].push_back(derivs.data() + n_qpoint*i_elem);
    }
  }
}

void Grid::populate_slice(std::vector<double>& elem_pos, std::vector<int> indices, int i_elem)
{
  if ((int)indices.size() < n_dim)
  {
    indices.push_back(0);
    for (int i = 0; i < basis.rank; ++i)
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
      stride /= basis.rank;
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

std::vector<double> Grid::get_pos(int i_elem)
{
  std::vector<double> elem_pos (n_qpoint*n_dim, 0.);
  std::vector<int> indices;
  populate_slice(elem_pos, indices, i_elem);
  return elem_pos;
}

double Grid::jacobian_det(int i_elem, int i_qpoint)
{
  return 1.;
}

bool Grid::execute_runge_kutta_stage()
{
  double* read0 = state_storage[0].data();
  double* read1 = state_w();
  double* write = (i_rk_stage == 2) ? read0 : read1;
  double weight1 = rk_weights[i_rk_stage]; double weight0 = 1. - weight1;
  for (int i = 0; i < n_elem*n_dof; ++i)
  {
    write[i] = weight1*read1[i] + weight0*read0[i];
  }
  ++i_rk_stage; ++i_write;
  if (i_write == 3) i_write = 1;
  if (i_rk_stage == 3) { i_rk_stage = 0; i_write = 1; ++iter; }
  i_read = i_rk_stage;
  return i_rk_stage == 0;
}

double Grid::get_stable_cfl()
{
  if ((basis.rank > 0) && (basis.rank <= 9))
  {
    return stable_cfl[basis.rank - 1];
  }
  else
  {
    throw std::runtime_error("Stable CFL number unknown for basis of desired rank.");
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
  for (int stride = n_qpoint/basis.rank, n_rows = 1; n_rows < n_qpoint;
       stride /= basis.rank, n_rows *= basis.rank)
  {
    for (int i_outer = 0; i_outer < n_rows; ++i_outer)
    {
      for (int i_inner = 0; i_inner < stride; ++i_inner)
      {
        for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
        {
          weights((i_outer*basis.rank + i_qpoint)*stride + i_inner)
          *= weights_1d(i_qpoint)*mesh_size;
        }
      }
    }
  }

  std::vector<double> total;
  double* sr = state_r();
  for (int i_elem = 0; i_elem < n_elem; ++i_elem)
  {
    std::vector<double> elem_pos = get_pos(i_elem);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
    {
      std::vector<double> point_pos;
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        point_pos.push_back(elem_pos[i_qpoint + i_axis*n_qpoint]);
      }
      std::vector<double> point_state;
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        point_state.push_back(sr[i_qpoint + i_var*n_qpoint + i_elem*n_dof]);
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
