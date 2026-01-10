#ifndef HEXED_ELEMENT_HPP_
#define HEXED_ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Kernel_element.hpp"
#include "Storage_params.hpp"
#include "Basis.hpp"
#include "Lock.hpp"
#include "reciprocal.hpp"
#include "Block.hpp"
#include "Face.hpp"

namespace hexed {

class Tree;
class Accessible_mesh;

/*! \brief Stores data associated with one mesh element.
 * \details Container only---does not have implementations of or information about the basis and algorithms.
 * This class represents a Cartesian (i.e., regular) element.
 * See also derived class `Deformed_element`.
 */
class Element : public Kernel_element, public Mortal {
  protected:
  // constructor that allows the vertices to be created as mobile, for the  benefit of `Deformed_element`
  Element(Storage_params, Tree&, bool mobile_vertices, int aniso_r_level, bool is_def);
  Mat<3> _compute_pos() const;
  Storage_params params;
  int n_dim;
  int _aniso_r_level;
  std::unique_ptr<next::Element_shape> _shape;

  private:
  int n_dof;
  int n_vert;
  int data_size;
  int face_size;
  Eigen::VectorXd data;
  Array<double> _vertex_data;
  std::array<double*, 6> faces; //!< layout: [2*i_dim + face_sign][i_var][i_qpoint]
  int _mask;
  // may contain a fake element that `this` is a subset of
  std::shared_ptr<next::Element_shape> _fake_shape;
  std::vector<Face> _faces;
  Array<int> _refinement_data;
  void _set_glued_pos();
  int _get_i_bf();
  friend Accessible_mesh; // necessary for `Accessible_mesh::set_mask`... need a better way to do this

  public:
  std::array<int, 6> face_record; //!< \brief for algorithms to book-keep information related to faces
  //! \brief Pointer to state data at faces. Must be populated by user
  double uncertainty = 0; //!< \brief refinement algorithms should set this value to some uncertainty metric
  static constexpr bool is_deformed = false; //!< \brief is this `Element` subclass deformed?
  Reciprocal_ptr<Element, Tree> tree; //!< \brief `Tree` this element was created from
  bool unrefinement_locked = false; //!< \brief if this is set to `true`, `Mesh_interface::update()` won't unrefine it
  //! \brief if `true`, this element has a face on the surface which was not properly snapped
  bool snapping_problem = false;
  //! \brief once any faces of this element have been snapped to the surface, set this to `false`
  bool needs_snapping = true;
  const Mat<> origin; //!< \brief origin which integer coordinates are relative to
  double residual;
  double flux_uncert;
  Lock lock; //!< \brief for any tasks where multiple threads might access an element simultaneously
  bool has_shock;
  bool spread_shock;

  Element(Storage_params, Tree& tree, int aniso_ref_level = 0);
  //! \details Can't copy an Element. Doing so would have to either duplicate or break vertex connections,
  //! both of which seem error prone.
  Element(const Element&) = delete;
  Element& operator=(const Element&) = delete;
  ~Element() = default;

  virtual inline bool get_is_deformed() {return is_deformed;} //!< for determining whether a pointer is deformed
  bool is_extruded();
  Storage_params storage_params() const;
  Array<double> position(const Basis&) const;
  Array<double> face_position(const Basis&) const;
  virtual void set_jacobian(const Basis& basis);
  double nominal_size() const;
  double nominal_shape(int i_dim) const override;
  double nominal_volume() const;
  int refinement_level();
  int aniso_ref_level();
  int& desired_refinement(int i_dim);
  int desired_refinement(int i_dim) const;
  Array<int> refinement_floor();
  Array<Int> nominal_position();
  double wall_distance() const; //!< \brief The distance from the farthest vertex to the wall.
  int wall_dimension();
  bool has_wall();
  bool is_sharp(int i_dim);
  //! pointer to state data for `i_stage`th Runge-Kutta stage.
  double* stage(int i_stage); //!< layout: [i_var][i_qpoint]
  Array<double> flow_state(); //!< layout: [i_var][i_qpoint]
  Array<double> numeric_state(); //!< layout: [i_var][i_qpoint]
  double* advection_state(); //!< layout: [i_node][i_qpoint] \note `0 <= i_node < row_size`
  //! pointer to scaling factor for local time step.
  double* time_step_scale() override; //!< layout: [i_qpoint]
  double* bulk_av_coef(); //!< layout: [i_qpoint]
  double* laplacian_av_coef(); //!< layout: [i_qpoint]
  double* art_visc_forcing(); //!< layout: [i_forcing][i_qpoint]
  //! \brief returns whether the element is included in the masked mesh.
  //! \details value can be set with `Accessible_mesh::set_mask`
  int mask() const override {return _mask;}
  Array<double> spectral_uncert(); //!< layout: [i_dim]

  /*! \brief Compute the Jacobian matrix.
   * \details I.e., derivative of `i_dim`th
   * physical coordinate wrt `j_dim`th reference coordinate. Trivial for this
   * class, may be non-trivial for derived (see `Deformed_element`).
   * For convenience, not performance.
   * For high-performance, use `double* Deformed_element::jacobian()`
   */
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint) const;
  virtual inline double jacobian_determinant(int i_qpoint) const {return 1.;} //!< \brief determinant of `jacobian`
  //! \brief Effective element dimensions accounting for deformation.
  //! \details layout: [i_dim]
  //! `mean_shape[i_dim]` is the mean over the element of the norm of the `i_dim`th column of the Jacobian matrix.
  Array<double> mean_shape(const Basis&) const;

  //! \brief Time step scale at the vertices. TSS in the interior is set by interpolating this.
  double& vertex_time_step_scale(int i_vertex) override;
  double& vertex_elwise_av(int i_vertex);
  double& vertex_fix_admis_coef(int i_vertex);
  void set_face(int i_face, double* data);
  bool is_connected(int i_face);
  inline Face& face(int i_face) {return _faces[i_face];}

  void create_shape(next::Mesh_blocks&, int boundary_face = next::Mesh_blocks::no_face);
  void create_fake(next::Mesh_blocks&);
  bool shared_fake() const;
  void split_shape(Element& split_from, double at, int from_face);
  void glue_shape(Element& glue_to, std::array<std::vector<double>, 2> corners);
  void destroy_shape();
  void destroy_fake();
  next::Element_shape& shape();
  const next::Element_shape& shape() const;
  // note `_fake_shape.use_count() == 0` does not imply `_fake_shape.get() == nullptr`
  inline next::Element_shape* fake_shape() {return _fake_shape.use_count() ? _fake_shape.get() : nullptr;}
  inline bool has_shape() const {return bool(_shape);}
  next::Element_shape& active_shape();

  double* state() override;
  double* residual_cache() override;
  double* face(int i_face, bool is_ldg) override;
  bool deformed() const override;
  double* reference_level_normals() override;
  double* jacobian_determinant() override;
  double* kernel_face_normal(int i_face) override;
  double* debug_variables() override;
  double& uncert() override;
};

}
#endif
