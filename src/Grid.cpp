#include <iostream>

#include <Grid.hpp>

Grid::Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg)
: n_var(n_var_arg), n_dim(n_dim_arg), n_elem(n_elem_arg), mesh_size(mesh_size_arg),
basis(basis_arg), iter(0), time(0.), i_rk_stage(0), i_read(0), i_write(1)
{
  n_qpoint = 1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    n_qpoint *= basis.rank;
  }
  n_dof = n_qpoint*n_var;
  for (int i = 0; i < (int)state_storage.size(); ++i) state_storage[i].resize(n_dof*n_elem, 0.);
  pos.resize(n_elem*n_dim, 0);
  for (int i_dim = 0; i_dim < 3*n_dim; ++i_dim)
  {
    std::vector<double*> empty;
    neighbor_storage.push_back(empty);
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
  std::array<double*, 3> state_data {state_storage[0].data(), state_storage[1].data(),
                                     state_storage[1].data()};
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
          if (pos_diff[i_dim] == 1)
          {
            for (int i_stage = 0; i_stage < 3; ++i_stage)
            {
              neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data[i_stage]
                                                                + n_dof*i_elem);
              neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data[i_stage]
                                                                + n_dof*j_elem);
            }
          }
          else if (pos_diff[i_dim] == -1)
          {
            for (int i_stage = 0; i_stage < 3; ++i_stage)
            {
              neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data[i_stage]
                                                                + n_dof*j_elem);
              neighbor_storage[i_dim + i_stage*n_dim].push_back(state_data[i_stage]
                                                                + n_dof*i_elem);
            }
          }
        }
      }
    }
  }
}

void Grid::clear_neighbors()
{
  for (int i_dim = 0; i_dim < 2*n_dim; ++i_dim)
  {
    neighbor_storage[i_dim].clear();
  }
}

void Grid::auto_connect()
{
  std::vector<int> periods;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) periods.push_back(0);
  auto_connect(periods);
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
                                           + pos[i_elem*n_dim + i_dim])*mesh_size;
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
  ++i_rk_stage;
  i_write = i_rk_stage + 1;
  if (i_write == 3) i_write = 1;
  if (i_rk_stage == 3) { i_rk_stage = 0; ++iter; }
  i_read = i_rk_stage;
  return i_rk_stage == 0;
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

Eigen::VectorXd Grid::state_integral()
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
  Eigen::Map<Eigen::MatrixXd> elem_mat (state_r(), n_qpoint, n_elem*n_var);
  Eigen::MatrixXd elem_integral = elem_mat.transpose()*weights;
  elem_integral.resize(n_var, n_elem);
  return elem_integral.rowwise().sum();
}
