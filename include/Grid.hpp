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
  std::vector<int> pos;
  double mesh_size;
  Basis& basis;
  int iter;
  double time;
  std::vector<double> origin;
  std::vector<Ghost_boundary_condition*> ghost_bound_conds;

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);

  // access
  virtual Element& element(int i_elem) = 0;
  virtual std::vector<double> get_pos(int i_elem) = 0;
  virtual double stable_time_step(double cfl_by_stable_cfl, Kernel_settings& setttings) = 0;
  int i_stage_read();
  int i_stage_write();

  // time integration
  bool execute_runge_kutta_stage();
  double get_stable_cfl();
  virtual void execute_write_face(Kernel_settings&) = 0;
  virtual void execute_neighbor(Kernel_settings&) = 0;
  virtual void execute_local(Kernel_settings&) = 0;
  virtual void execute_req_visc(Kernel_settings&) = 0;
  virtual void execute_cont_visc(Kernel_settings&) = 0;
  virtual void execute_local_derivative(int i_var, int i_dim, Kernel_settings&) = 0;
  virtual void execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings&) = 0;
  virtual void execute_av_flux(Kernel_settings&) = 0;
  virtual void execute_local_av(int i_var, int i_dim, Kernel_settings&) = 0;
  virtual void execute_neighbor_av(int i_var, int i_dim, Kernel_settings&) = 0;

  // modification
  virtual int add_element(std::vector<int> position);

  // diagnostic
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
};

}
#endif
