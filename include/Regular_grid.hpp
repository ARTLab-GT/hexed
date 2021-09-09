#ifndef CARTDG_REGULAR_GRID_HPP_
#define CARTDG_REGULAR_GRID_HPP_

#include "Grid.hpp"

namespace cartdg
{

class Regular_grid : public Grid
{
  public:
  Regular_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);

  virtual Element& element(int i_elem);
  elem_con connection(int i_dim, int i_con); // mostly for testing
  int n_con(int i_dim); // mostly for testing
  virtual std::vector<double> get_pos(int i_elem);

  virtual double jacobian_det(int i_elem, int i_qpoint);
  std::vector<double**> neighbor_connections_r(); // FIXME: return a pointer
  std::vector<double**> neighbor_connections_w();
  std::vector<double**> deriv_neighbor_connections();
  std::vector<double**> visc_neighbor_connections();
  std::vector<int> n_neighb_con();
  virtual double stable_time_step(double cfl_by_stable_cfl, Kernel_settings& setttings);

  virtual void execute_local(Kernel_settings&);
  virtual void execute_neighbor(Kernel_settings&);
  virtual void execute_req_visc(Kernel_settings&);
  virtual void execute_cont_visc(Kernel_settings&);
  virtual void execute_local_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_av_flux(Kernel_settings&);
  virtual void execute_local_av(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_av(int i_var, int i_dim, Kernel_settings&);

  virtual int add_element(std::vector<int> position);
  void add_connection(int i_elem0, int i_elem1, int i_dim);
  virtual void clear_neighbors();
  void auto_connect(std::vector<int> periods);
  void auto_connect();

  private:
  elem_vec elements;
  elem_con_vec elem_cons;
  void populate_slice(std::vector<double>&, std::vector<int>, int);
};

}
#endif
