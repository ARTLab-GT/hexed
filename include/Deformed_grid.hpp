#ifndef CARTDG_DEFORMED_GRID_HPP_
#define CARTDG_DEFORMED_GRID_HPP_

#include "Grid.hpp"
#include "Vertex.hpp"

namespace cartdg
{

class Deformed_grid : public Grid
{
  public:
  std::vector<Vertex> vertices;
  std::vector<int> vertex_ids;
  std::vector<double> node_adjustments;
  std::vector<double> jacobian;
  std::vector<double> default_jacobian;
  std::vector<int> neighbor_inds;
  std::vector<double*> jacobian_neighbors;
  std::vector<int> neighbor_axes;
  std::vector<int> neighbor_is_positive;
  int n_vertices;
  std::vector<int> i_elem_wall;
  std::vector<int> i_dim_wall;
  std::vector<int> is_positive_wall;

  Deformed_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg,
                Basis& basis_arg);
  Vertex& get_vertex(int i_vertex);
  virtual int add_element(std::vector<int> position);
  virtual std::vector<double> get_pos(int i_elem);
  void add_wall(int i_elem, int i_dim, bool is_positive_face);
  virtual double jacobian_det(int i_elem, int i_qpoint);
  virtual void execute_local(Kernel_settings&);
  virtual void execute_neighbor(Kernel_settings&);
  virtual void execute_req_visc(Kernel_settings&);
  virtual void execute_cont_visc(Kernel_settings&);
  virtual void execute_local_derivative(int i_var, int i_dim, Kernel_settings&); // FIXME: for now, does nothing
  virtual void execute_neighbor_derivative(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_av_flux(Kernel_settings&);
  virtual void execute_local_av(int i_var, int i_dim, Kernel_settings&);
  virtual void execute_neighbor_av(int i_var, int i_dim, Kernel_settings&);

  // Note: the following functions must be called in the order
  // that they appear.
  void connect(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
               std::array<bool, 2> is_positive);
  void calc_jacobian();
  void update_connections();
  void connect_non_def(std::array<int, 2> i_elem, std::array<int, 2> i_dim,
                       std::array<bool, 2> is_positive, Grid& other_grid);

  virtual void visualize(std::string file_name);
  std::vector<double> face_integral(Domain_func& integrand, int i_elem, int i_dim, bool is_positive);
  std::vector<double> surface_integral(Domain_func& integrand);
  double** state_connections_r();
  double** state_connections_w();

  private:
  void add_vertices(std::vector<int> position, int i_dim);
};

}

#endif
