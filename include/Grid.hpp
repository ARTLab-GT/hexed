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
  Basis& basis;
  int iter;
  double time;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
  virtual ~Grid();

  double* state_r();
  double* state_w();
  std::vector<double**> neighbor_connections_r();
  std::vector<double**> neighbor_connections_w();
  std::vector<int> n_neighb_con();
  std::vector<double> get_pos(int i_elem);
  bool execute_runge_kutta_stage();

  void auto_connect(std::vector<int> periods);
  void auto_connect();
  void clear_neighbors();
  void visualize(std::string file_name);
  void print();
  Eigen::VectorXd state_integral();

  protected:
  int i_rk_stage;
  int i_read;
  int i_write;
  double rk_weights [3] {1., 1./4., 2./3.};
  std::array<std::vector<double>, 3> state_storage {};
  std::vector<std::vector<double*>> neighbor_storage;

  private:
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

#endif
