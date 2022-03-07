#ifndef CARTDG_GRID_HPP_
#define CARTDG_GRID_HPP_

#include <vector>
#include <string>

#include <Eigen/Dense>

#include "Basis.hpp"
#include "Element.hpp"
#include "Kernel_settings.hpp"
#include "Domain_func.hpp"
#include "Ghost_boundary_condition.hpp"
#include "Hanging_vertex_matcher.hpp"

namespace cartdg
{

class Tecplot_file;

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

  Grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, Basis& basis_arg);

  // access
  virtual Element& element(int i_elem) = 0;
  virtual std::vector<double> get_pos(int i_elem) = 0;
  // maximum characteristic speed in reference space
  virtual double max_reference_speed(Kernel_settings& setttings) = 0;
  int i_stage_read();
  int i_stage_write();

  // time integration
  bool execute_runge_kutta_stage();
  virtual void execute_write_face(Kernel_settings&) = 0;
  virtual void execute_neighbor(Kernel_settings&) = 0;
  virtual void execute_local(Kernel_settings&) = 0;
  virtual double execute_req_visc(Kernel_settings&) = 0;
  virtual void execute_write_face_gradient(int i_var, Kernel_settings&) = 0;
  virtual void execute_neighbor_gradient(int i_var, Kernel_settings&) = 0;
  virtual void execute_local_gradient(int i_var, Kernel_settings&) = 0;
  virtual void execute_write_face_av(int i_var, Kernel_settings&) = 0;
  virtual void execute_neighbor_av(int i_var, Kernel_settings&) = 0;
  virtual void execute_local_av(int i_var, Kernel_settings&) = 0;

  // modification
  virtual int add_element(std::vector<int> position);
  virtual void add_element_gbc(int i_elem, Ghost_boundary_condition&) = 0;
  void match_hanging(Element::shareable_value_access);

  // diagnostic
  void print();

  protected:
  int i_rk_stage;
  int i_read;
  int i_write;
  double rk_weights [3] {1., 1./4., 2./3.};
  Storage_params storage_params;
  std::vector<Hanging_vertex_matcher> hanging_matchers;
};

}
#endif
