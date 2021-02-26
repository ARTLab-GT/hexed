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
  std::vector<int> neighbor_inds;
  std::vector<double*> jacobian_neighbors;
  std::vector<int> neighbor_axes;
  std::vector<int> neighbor_is_positive;
  std::vector<int> neighbor_is_deformed;
  int n_vertices;

  Deformed_grid(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg,
                Basis& basis_arg);
  Vertex& get_vertex(int i_vertex);
  virtual int add_element(std::vector<int> position);
  virtual std::vector<double> get_pos(int i_elem);
  void connect(std::array<int, 2> i_elem, std::array<int, 2> i_axis,
               std::array<bool, 2> is_positive);
  void update_connections();
  virtual void visualize(std::string file_name);
  void calc_jacobian();
  double** state_connections_r();
  double** state_connections_w();

  private:
  void add_vertices(std::vector<int> position, int i_dim);
};

}

#endif