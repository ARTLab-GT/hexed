#ifndef GRID_HPP_
#define GRID_HPP_

#include <vector>
#include <string>

#include <Eigen/Dense>

#include <Basis.hpp>

class Grid
{
  public:
  int n_var;
  int n_dim;
  int n_qpoint;
  int n_dof;
  int n_elem;
  std::vector<int> pos;
  double mesh_size;
  double time;
  Basis& basis;
  int iter;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
  virtual ~Grid();

  inline double* state_r() { return ((iter%2 == 0) ? state_r_storage : state_w_storage).data(); }
  inline double* state_w() { return ((iter%2 == 1) ? state_r_storage : state_w_storage).data(); }
  std::vector<double**> neighbor_connections();
  std::vector<int> n_neighb_con();
  std::vector<double> get_pos(int i_elem);

  void auto_connect(std::vector<int> periods);
  void auto_connect();
  void clear_neighbors();
  void visualize(std::string file_name);
  void print();
  Eigen::VectorXd state_integral();

  protected:
  std::vector<double> state_r_storage;
  std::vector<double> state_w_storage;
  std::vector<std::vector<double*>> neighbor_storage;

  private:
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

#endif
