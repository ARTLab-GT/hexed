#ifndef DEFORMED_ELEMENT_HPP_
#define DEFORMED_ELEMENT_HPP_

#include "Element.hpp"
#include "Vertex.hpp"

namespace cartdg
{

/*
 * Represents an Element which is not a perfect square/cube. Note: Jacobian matrix
 * is nontrivial.
 */
class Deformed_element : public Element
{
  int n_qpoint;
  Eigen::VectorXd jac;
  std::vector<Vertex::Transferable_ptr> vertices;
  Eigen::VectorXd node_adj;

  public:
  bool degenerate = 0;

  /*
   * `pos` sets the position of vertex 0 in increments of `mesh_size`. Unspecified
   * components set to 0. `mesh_size` sets distance between vertices. Vertex positions
   * may be modified later.
   */
  Deformed_element(Storage_params, std::vector<int> pos={}, double mesh_size=1.);
  // Pointer to jacobian data. Use this for performance-cricial applications.
  double* jacobian(); // Layout: [i_dim][j_dim][i_qpoint]
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  inline Vertex& vertex(int i_vertex) { return *vertices[i_vertex]; }
  // No simple way to explain what this represents.
  double* node_adjustments(); // Layout: [i_dim][is_positive][i_face_qpoint]
};

class Deformed_elem_con
{
  public:
  std::array<Deformed_element*, 2> element;
  // Following two members designate which faces are participating in connection
  std::array<int, 2> i_dim;
  std::array<bool, 2> is_positive;
};

// Designates that a face of an element is participating in a wall boundary condition
class Deformed_elem_wall
{
  public:
  Deformed_element* element;
  int i_dim;
  bool is_positive;
  int i_elem;
};

typedef std::vector<std::unique_ptr<Deformed_element>> def_elem_vec;
typedef std::vector<Deformed_elem_con> def_elem_con_vec;
typedef std::vector<Deformed_elem_wall> def_elem_wall_vec;

}
#endif
