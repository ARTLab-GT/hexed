#ifndef CARTDG_GRID_HPP_
#define CARTDG_GRID_HPP_

#include <vector>
#include <string>

#include <Eigen/Dense>

#include "Basis.hpp"
#include "Element.hpp"
#include "Ghost_boundary_condition.hpp"
#include "Kernel_settings.hpp"
#include "Domain_func.hpp"

namespace cartdg
{

class Grid
{
  public:
  int n_var;
  int n_dim;
  int n_vertices;
  int n_qpoint;
  int n_dof;
  int n_elem;
  std::array<std::vector<double>, 3> state_storage {};
  std::vector<std::vector<double*>> neighbor_storage;
  std::vector<double> derivs;
  std::vector<std::vector<double*>> deriv_neighbor_storage;
  std::vector<double> visc;
  std::vector<std::vector<double*>> visc_neighbor_storage;
  std::vector<int> viscous_inds;
  std::vector<int> pos;
  double mesh_size;
  Basis& basis;
  int iter;
  double time;
  std::vector<double> origin;
  std::vector<Ghost_boundary_condition*> ghost_bound_conds;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);
  Grid(Grid&&) = default;
  virtual ~Grid();

  // functions for accessing data
  Element& element(int i_elem);
  double* state_r();
  double* state_w();
  std::vector<double**> neighbor_connections_r(); // FIXME: return a pointer
  std::vector<double**> neighbor_connections_w();
  std::vector<double**> deriv_neighbor_connections();
  std::vector<double**> visc_neighbor_connections();
  std::vector<int> n_neighb_con();
  virtual std::vector<double> get_pos(int i_elem);
  virtual double jacobian_det(int i_elem, int i_qpoint);
  virtual double stable_time_step(double cfl_by_stable_cfl, Kernel_settings& setttings);

  // functions that execute some aspect of time integration
  bool execute_runge_kutta_stage();
  double get_stable_cfl();
  virtual void execute_local(Kernel_settings&);
  virtual void execute_neighbor(Kernel_settings&);
  virtual void execute_req_visc(Kernel_settings&);
  virtual void execute_cont_visc(Kernel_settings&);
  virtual void execute_local_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_av_flux(Kernel_settings&);
  virtual void execute_local_av(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_av(int i_var, int i_dim, Kernel_settings&);

  // functions that resize/reallocate/modify data
  void auto_connect(std::vector<int> periods);
  void auto_connect();
  void clear_neighbors();
  virtual int add_element(std::vector<int> position);
  void add_connection(int i_elem0, int i_elem1, int i_dim);

  // functions that provide diagnostic information
  virtual void visualize(std::string file_name);
  void print();
  std::vector<double> integral();
  std::vector<double> integral(Domain_func& integrand);

  protected:
  int i_rk_stage;
  int i_read;
  int i_write;
  double rk_weights [3] {1., 1./4., 2./3.};
  double stable_cfl [9] {1.256, 0.409, 0.209, 0.130, 0.089, 0.066, 0.051, 0.040, 0.033};
  Storage_params storage_params;
  elem_vec elements;

  private:
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

}
#endif
