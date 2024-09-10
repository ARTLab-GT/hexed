#ifndef HEXED_ELEMENT_HPP_
#define HEXED_ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Kernel_element.hpp"
#include "Storage_params.hpp"
#include "Vertex.hpp"
#include "Basis.hpp"
#include "Lock.hpp"
#include "Mutual_ptr.hpp"
#include "Block.hpp"

namespace hexed {

class Tree;
class Accessible_mesh;

/*! \brief Stores data associated with one mesh element.
 * \details Container only---does not have implementations of or information about the basis and algorithms.
 * This class represents a Cartesian (i.e., regular) element.
 * See also derived class `Deformed_element`.
 */
class Element : public Kernel_element {
  protected:
  Storage_params params;
  int n_dim;
  std::vector<int> _nom_pos;
  double _nom_sz;
  int _r_level;
  int _aniso_r_level;
  std::vector<Vertex::Transferable_ptr> vertices;
  std::unique_ptr<next::Element_shape> _shape;
  // constructor that allows the vertices to be created as mobile, for the  benefit of `Deformed_element`
  Element(Storage_params, std::vector<int> pos, double mesh_size, int ref_level, Mat<> origin_arg, bool mobile_vertices, int aniso_r_level);

  private:
  int n_dof;
  int n_vert;
  int data_size;
  Eigen::VectorXd data;
  Eigen::VectorXd vertex_data;
  std::array<double*, 6> faces; //!< layout: [2*i_dim + face_sign][i_var][i_qpoint]
  int _mask;
  friend Accessible_mesh; // necessary for `Accessible_mesh::set_mask`... need a better way to do this

  public:
  std::array<int, 6> face_record; //!< \brief for algorithms to book-keep information related to faces
  //! \brief Pointer to state data at faces. Must be populated by user
  double uncertainty = 0; //!< \brief refinement algorithms should set this value to some uncertainty metric
  static constexpr bool is_deformed = false; //!< \brief is this `Element` subclass deformed?
  Mutual_ptr<Element, Tree> tree; //!< \brief `Tree` this element was created from
  bool unrefinement_locked = false; //!< \brief if this is set to `true`, `Mesh_interface::update()` won't unrefine it
  bool snapping_problem = false; //!< \brief if `true`, this element has a face on the surface which was not properly snapped
  bool needs_snapping = true; //!< \brief once any faces of this element have been snapped to the surface, set this to `false`
  const Mat<> origin; //!< \brief origin which integer coordinates are relative to
  Lock lock; //!< \brief for any tasks where multiple threads might access an element simultaneously

  /*!
   * The `Storage_params` defines the amount of storage that must be allocated.
   * `pos` specifies the position of vertex 0 relative to `origin_arg` in intervals of the nominal size.
   * The nominal size is defined to be `mesh_size`/(2^`ref_level`).
   * The vertices will be spaced at intervals of the nominal size.
   * Only the first `n_dim` elements of `origin_arg` are considered.
   */
  Element(Storage_params, std::vector<int> pos = {}, double mesh_size = 1., int ref_level = 0,
          Mat<> origin_arg = Mat<>::Zero(3), int aniso_ref_level = 0);
  virtual inline bool get_is_deformed() {return is_deformed;} //!< for determining whether a pointer is deformed
  //! Can't copy an Element. Doing so would have to either duplicate or break vertex connections, both of which seem error prone.
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;
  ~Element() = default;

  Storage_params storage_params();
  Array<double> position(const Basis&) const;
  Array<double> vis_position(const Basis&) const;
  Array<double> face_position(const Basis&) const;
  virtual void set_jacobian(const Basis& basis);
  inline double nominal_size() const override {return _nom_sz;}
  inline int refinement_level() {return _r_level;} //!< \brief indicates how many times this element has been isotropically refined
  int aniso_ref_level() {return _aniso_r_level;} //!< \brief indicates how many times this element has been anisotropically refined
  inline std::vector<int> nominal_position() {return _nom_pos;}
  //! pointer to state data for `i_stage`th Runge-Kutta stage.
  double* stage(int i_stage); //!< layout: [i_var][i_qpoint]
  double* advection_state(); //!< layout: [i_node][i_qpoint] \note `0 <= i_node < row_size`
  //! pointer to scaling factor for local time step.
  double* time_step_scale() override; //!< layout: [i_qpoint]
  double* bulk_av_coef(); //!< layout: [i_qpoint]
  double* laplacian_av_coef(); //!< layout: [i_qpoint]
  double* art_visc_forcing(); //!< layout: [i_forcing][i_qpoint]
  //! \brief returns whether the element is included in the masked mesh.
  //! \details value can be set with `Accessible_mesh::set_mask`
  int mask() const override {return _mask;}

  /*! \brief Compute the Jacobian matrix.
   * \details I.e., derivative of `i_dim`th
   * physical coordinate wrt `j_dim`th reference coordinate. Trivial for this
   * class, may be non-trivial for derived (see `Deformed_element`).
   * For convenience, not performance.
   * For high-performance, use `double* Deformed_element::jacobian()`
   */
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  virtual inline double jacobian_determinant(int i_qpoint) {return 1.;} //!< \brief determinant of `jacobian`

  inline Vertex& vertex(int i_vertex) { return *vertices[i_vertex]; }
  template <int i_dim> double& vertex_position(int i_vertex) {return vertices[i_vertex]->pos[i_dim];}
  //! \brief functions to communicate with nodal neighbors
  void push_shareable_value(std::function<double(Element&, int i_vertex)>); //!< \brief writes shareable value to vertices so that shared value can be determined
  void fetch_shareable_value(std::function<double&(Element&, int i_vertex)> access_fun, std::function<double(Mat<>)> reduction = Vertex::vector_max); // set `this`'s copy of shareable value to the shared values at the vertices
  //! \brief Time step scale at the vertices. TSS in the interior is set by interpolating this.
  double& vertex_time_step_scale(int i_vertex) override;
  double& vertex_elwise_av(int i_vertex);
  double& vertex_fix_admis_coef(int i_vertex);
  void set_needs_smooth(bool); //!< \brief sets the `Vertex::Transferable_ptr::needs_smooth` of the vertices
  void set_face(int i_face, double* data);
  bool is_connected(int i_face);

  void create_shape(next::Mesh_blocks&, int boundary_face = next::Mesh_blocks::no_face);
  next::Element_shape& shape();

  double* state() override;
  double* residual_cache() override;
  double* face(int i_face, bool is_ldg) override;
  bool deformed() const override;
  double* reference_level_normals() override;
  double* jacobian_determinant() override;
  double* kernel_face_normal(int i_face) override;
  double& uncert() override;
};

}
#endif
