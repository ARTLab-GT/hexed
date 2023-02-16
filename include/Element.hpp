#ifndef HEXED_ELEMENT_HPP_
#define HEXED_ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Storage_params.hpp"
#include "Vertex.hpp"
#include "Basis.hpp"

namespace hexed
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
  double nom_sz;
  int r_level;
  std::vector<Vertex::Transferable_ptr> vertices;
  void set_wall_dist(const Basis& basis);

  private:
  int n_dof;
  int n_vert;
  int data_size;
  Eigen::VectorXd data;
  Eigen::VectorXd vertex_tss;

  public:
  // pointer to a function that can access some data associated with the vertices
  typedef double& (Element::*vertex_value_access)(int i_vertex);
  std::array<int, 6> face_record; // for algorithms to book-keep information related to faces
  std::array<double*, 6> faces; // layout: [2*i_dim + face_sign][i_var][i_qpoint]
  double resolution_badness = 0;
  static constexpr bool is_deformed = false;

  /*
   * The `Storage_params1 defines the amount of storage that must be allocated.
   * `pos` specifies the position of vertex 0 in intervals of the nominal size.
   * The nominal size is defined to be `mesh_size`/(2^`ref_level`).
   * The vertices will be spaced at intervals of the nominal size.
   */
  Element(Storage_params, std::vector<int> pos={}, double mesh_size=1., int ref_level = 0);
  virtual inline bool get_is_deformed() {return is_deformed;} // for determining whether a pointer is deformed
  // Can't copy an Element. Doing so would have to either duplicate or break vertex connections,
  // both of which seem error prone.
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;
  ~Element() = default;

  Storage_params storage_params();
  virtual std::vector<double> position(const Basis&, int i_qpoint); // note: ignores vertex positions
  // obtains face position based on interior qpoint positions (as defined by `position()`)
  std::vector<double> face_position(const Basis&, int i_face, int i_face_qpoint);
  virtual void set_jacobian(const Basis& basis);
  inline double nominal_size() {return nom_sz;}
  inline int refinement_level() {return r_level;}
  inline std::vector<int> nominal_position() {return nom_pos;}
  // Pointer to state data for `i_stage`th Runge-Kutta stage.
  double* stage(int i_stage); // Layout: [i_var][i_qpoint]
  // pointer to advection state data
  double* advection_state();
  // pointer to scaling factor for local time step.
  double* time_step_scale(); // Layout: [i_qpoint]
  double* art_visc_coef(); // layout: [i_qpoint]
  double* art_visc_forcing(); // layout: [i_qpoint]
  // Pointer state data at faces. Must be populated by user
  virtual double* node_adjustments() {return nullptr;} // overriden by `Deformed_element`

  /*
   * Following two functions compute the Jacobian. I.e., derivative of `i_dim`th
   * physical coordinate wrt `j_dim`th reference coordinate. Trivial for this
   * class, may be non-trivial for derived (see `Deformed_element`). For
   * convenience, not performance.
   */
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint); // identity matrix
  virtual inline double jacobian_determinant(int i_qpoint) {return 1.;}

  inline Vertex& vertex(int i_vertex) { return *vertices[i_vertex]; }
  template <int i_dim> double& vertex_position(int i_vertex) {return vertices[i_vertex]->pos[i_dim];}
  // functions to communicate with nodal neighbors
  void push_shareable_value(vertex_value_access access_func); // writes shareable value to vertices so that shared value can be determined
  void fetch_shareable_value(vertex_value_access access_func, Vertex::reduction = Vertex::vector_max); // set `this`'s copy of shareable value to the shared values at the vertices
  // Time step scale at the vertices. TSS in the interior is set by interpolating this.
  double& vertex_time_step_scale(int i_vertex);
};

}
#endif
