#ifndef CARTDG_DEFORMED_GRID_HPP_
#define CARTDG_DEFORMED_GRID_HPP_

#include "Grid.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

class Deformed_grid : public Grid
{
  def_elem_vec elements;
  std::vector<Vertex::Non_transferable_ptr> vertices;
  def_elem_con_vec elem_cons;
  def_reg_con_vec def_reg_cons;
  def_elem_wall_vec walls;

  public:
  Deformed_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg,
                Basis& basis_arg);

  // access
  virtual Element& element(int i_elem);
  Deformed_element& deformed_element(int i_elem); // points to same object as element(int) but with different type
  virtual double stable_time_step(double cfl_by_stable_cfl, Kernel_settings& setttings);
  // the following are mostly for testing
  Deformed_elem_con connection(int i_con);
  def_reg_con def_reg_connection(int i_dim, int i_con);
  Deformed_elem_wall def_elem_wall(int i_wall);

  // modification
  virtual int add_element(std::vector<int> position);
  virtual std::vector<double> get_pos(int i_elem);
  void add_wall(int i_elem, int i_dim, bool is_positive_face);
  void connect(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
               std::array<bool, 2> is_positive);
  void connect_non_def(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
                       std::array<bool, 2> is_positive, Grid& other_grid);
  void purge_vertices();
  void calc_vertex_relaxation();
  void apply_vertex_relaxation();
  void calc_jacobian(); // must be called after vertex locations are correct

  // time integration
  virtual void execute_local(Kernel_settings&);
  virtual void execute_neighbor(Kernel_settings&);
  virtual void execute_req_visc(Kernel_settings&);
  virtual void execute_cont_visc(Kernel_settings&);
  virtual void execute_local_derivative(int i_var, int i_dim, Kernel_settings&); // FIXME: for now, does nothing
  virtual void execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_av_flux(Kernel_settings&);
  virtual void execute_local_av(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_av(int i_var, int i_dim, Kernel_settings&);

  // diagnostic
  virtual void visualize_qpoints(std::string file_name);
  std::vector<double> face_integral(Domain_func& integrand, int i_elem, int i_dim, bool is_positive);
  std::vector<double> surface_integral(Domain_func& integrand);
};

}

#endif
