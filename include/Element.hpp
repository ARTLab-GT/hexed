#ifndef CARTDG_ELEMENT_HPP_
#define CARTDG_ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Storage_params.hpp"
#include "Vertex.hpp"
#include "Basis.hpp"

namespace cartdg
{

/*
 * Stores data associated with one element. Container only --
 * does not have implementations of or information about the basis and algorithms.
 * This class represents a Cartesian (i.e., regular) element. See also derived class
 * `Deformed_element`.
 */
class Element
{
  protected:
  Storage_params params;
  int n_dim;
  std::vector<int> nom_pos;
  std::vector<Vertex::Transferable_ptr> vertices;

  private:
  int n_dof;
  int n_vert;
  int data_size;
  double nom_sz;
  Eigen::VectorXd data;
  Eigen::VectorXd visc_storage;
  Eigen::VectorXd vertex_tss;
  Eigen::VectorXd derivative_storage;

  public:
  // this type represents a pointer to a function that can access some data member which is recorded at the vertices
  typedef double* (Element::*shareable_value_access)();

  Element(Storage_params, std::vector<int> pos={}, double mesh_size=1.);
  // Can't copy an Element. Doing so would have to either duplicate or break vertex connections,
  // both of which seem error prone.
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;
  ~Element() = default;

  Storage_params storage_params();
  virtual std::vector<double> position(const Basis&, int i_qpoint); // note: ignores vertex positions
  inline double nominal_size() {return nom_sz;}
  // Pointer to state data for `i_stage`th Runge-Kutta stage.
  double* stage(int i_stage); // Layout: [i_var][i_qpoint]
  // pointer to scaling factor for local time step.
  double* time_step_scale(); // Layout: [i_qpoint]
  // Pointer state data at faces. Must be populated by user
  double* face(); // Layout: [i_dim][is_positive][i_var][i_face_qpoint]

  /*
   * Following two functions compute the Jacobian. I.e., derivative of `i_dim`th
   * physical coordinate wrt `j_dim`th reference coordinate. Trivial for this
   * class, may be non-trivial for derived (see `Deformed_element`). For
   * convenience, not performance.
   */
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint); // identity matrix
  double jacobian_determinant(int i_qpoint); // returns 1.

  inline Vertex& vertex(int i_vertex) { return *vertices[i_vertex]; }
  // functions to communicate with nodal neighbors
  void push_shareable_value(shareable_value_access access_func); // writes shareable value to vertices so that shared value can be determined
  void fetch_shareable_value(shareable_value_access access_func, Vertex::reduction = Vertex::vector_max); // set `this`'s copy of shareable value to the shared values at the vertices

  double* viscosity(); // Artificial viscosity coefficient at corners.
  bool viscous(); // Should artificial viscosity be applied in this element?
  // Time step scale at the vertices. TSS in the interior is set by interpolating this.
  double* vertex_time_step_scale();
  // Pointer to storage for derivative. Must be populated by user.
  double* derivative(); // Layout: [i_qpoint]
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::array<double*, 2> elem_con; // pointers to the faces which are being connected
typedef std::vector<std::vector<elem_con>> elem_con_vec;

}

#endif
