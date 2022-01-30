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

class Deformed_face
{
  int n_dim;
  int n_fqpoint; // number of quadrature points on face
  Eigen::VectorXd jac;

  public:
  Deformed_face(Storage_params);
  /*
   * Pointer to a Jacobian matrix at the face mutually agreed upon by any
   * elements involved. Use this for performance-critical access to jacobian.
   * Layout: [i_dim][j_dim][i_face_qpoint]
   */
  double* jacobian();
  /*
   * element of the Jacobian matrix at row `i_dim` and column `j_dim`,
   * i.e. derivative of `i_dim`th physical coordinate wrt `j_dim`th reference coordinate,
   * evaluated at the `i_qpoint`th face quadrature point. For convenience, not
   * performance.
   */
  double jacobian(int i_dim, int j_dim, int i_qpoint);
};

// Stores the information required to identify a face of a deformed element
class Face_index
{
  public:
  Deformed_element* element;
  int i_dim; // which reference coordinate axis is normal to face
  bool is_positive; // `true` if `i_dim`th reference coordinate is 1.0 at face (`false` if 0.0)
};

// represents a connection between two deformed elements
class Deformed_elem_con
{
  std::array<Face_index, 2> face_inds;

  public:
  Deformed_elem_con(std::array<Face_index, 2>);
  Face_index face_index(int i_side);
  /*
   * return `true` if the face normal direction defined by the face jacobian is the opposite
   * of that defined by the jacobian of element `i_side`. Valid args: 0, 1
   */
  bool flip(int i_side);
};

// Designates that a face of an element is participating in a wall boundary condition
class Deformed_elem_wall
{
  Face_index f_ind;
  int i_el;

  public:
  Deformed_elem_wall(Face_index, int i_elem_arg);
  Face_index face_index();
  int i_elem();
};

typedef std::vector<std::unique_ptr<Deformed_element>> def_elem_vec;
typedef std::vector<Deformed_elem_con> def_elem_con_vec;
typedef std::vector<Deformed_elem_wall> def_elem_wall_vec;

}
#endif
