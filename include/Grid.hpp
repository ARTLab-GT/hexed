#ifndef CARTDG_GRID_HPP_
#define CARTDG_GRID_HPP_

#include <vector>
#include <string>

#include <Eigen/Dense>

#include "Basis.hpp"
#include "Fitted_boundary_condition.hpp"
#include "kernels/kernel_types.hpp"

namespace cartdg
{

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
  std::vector<double> origin;
  std::vector<Fitted_boundary_condition*> fit_bound_conds;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
  virtual ~Grid();

  // functions for accessing data
  double* state_r();
  double* state_w();
  std::vector<double**> neighbor_connections_r();
  std::vector<double**> neighbor_connections_w();
  std::vector<int> n_neighb_con();
  std::vector<double> get_pos(int i_elem);

  // functions that execute some aspect of time integration
  bool execute_runge_kutta_stage();
  void apply_fit_bound_conds(double d_t_by_d_x);
  double get_stable_cfl();

  // functions that resize/reallocate/modify data
  void auto_connect(std::vector<int> periods);
  void auto_connect();
  void clear_neighbors();
  int add_element(std::vector<int> position);
  void add_connection(int i_elem0, int i_elem1, int i_dim);

  // functions that provide diagnostic information
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
  double stable_cfl [9] {1.256, 0.409, 0.209, 0.130, 0.089, 0.066, 0.051, 0.040, 0.033};

  virtual Flux_kernel get_flux_kernel();
  virtual Read_kernel get_read_kernel();
  virtual Write_kernel get_write_kernel();

  private:
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

}
#endif
