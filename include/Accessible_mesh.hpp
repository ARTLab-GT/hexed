#ifndef HEXED_ACCESSIBLE_MESH_HPP_
#define HEXED_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Mesh_by_type.hpp"
#include "Tree.hpp"

namespace hexed
{

/*!
 * A mesh that supports access to the actual elements with the numerical data they contain. This level of
 * access is required by the numerical scheme but should be hidden from the library user, who should not be
 * concerned with numerical details.
 */
class Accessible_mesh : public Mesh
{
  Storage_params params;
  int n_vert;
  double root_sz;
  Mesh_by_type<         Element> car;
  Mesh_by_type<Deformed_element> def;
  // create a `Vector_view` that can look at `def.elements()` as `Element&`s.
  Vector_view<Element&, Deformed_element&, &trivial_convert<Element&, Deformed_element&>, Sequence> def_as_car;
  Concatenation<Element&> elems;
  Concatenation<Element_connection&> elem_cons;
  std::vector<Boundary_condition> bound_conds;
  std::unique_ptr<Surface_geom> surf_geom;
  int surf_bc_sn;
  Concatenation<Face_connection<Deformed_element>&> bound_face_cons;
  Concatenation<Boundary_connection&> bound_cons;
  Concatenation<Face_connection<Deformed_element>&> def_face_cons;
  Concatenation<Refined_face&> ref_face_v;
  Concatenation<Hanging_vertex_matcher&> matcher_v;
  std::vector<Vertex::Non_transferable_ptr> vert_ptrs;
  std::vector<int> extrude_cons;
  std::unique_ptr<Tree> tree; // could be null! don't forget to check
  std::vector<int> tree_bcs;

  Element_container& container(bool is_deformed);
  Element& add_elem(bool is_deformed, Tree&);
  bool is_surface(Tree*);
  template<typename element_t> Mesh_by_type<element_t>& mbt(); // gets either `car` or `def`
  template<typename element_t> void connect_new(int start_at); // connects new elements in `mbt<element_t>()`. helper function for `refine`
  void refine_by_record(bool is_deformed, int start, int end);
  bool needs_refine(Tree*);
  void deform();
  void purge();

  public:
  /*!
   * \param params parameters specifying what data is stored in each Element (row size, number of dimensions, etc.)
   * \param root_size defines the \ref root_size of the mesh.
   */
  Accessible_mesh(Storage_params params, double root_size);
  inline double root_size() override {return root_sz;}
  //! \returns a View_by_type containing only the Cartesian elements in the mesh
  inline View_by_type<         Element>& cartesian() {return car;}
  //! \returns a View_by_type containing only the deformed elements in the mesh
  inline View_by_type<Deformed_element>&  deformed() {return def;}
  int add_element(int ref_level, bool is_deformed, std::vector<int> position) override;
  //! Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
  //! access all elements, both Cartesian and deformed
  Sequence<Element&>& elements() {return elems;}
  void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                         std::array<bool, 2> is_deformed = {false, false}) override;
  void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction) override;
  void connect_hanging(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Deformed_element>,
                               bool coarse_deformed = false, std::vector<bool> fine_deformed = {false, false, false, false},
                               std::array<bool, 2> stretch = {false, false}) override;
  //! \returns a view of all connections between elements, including one connection for every fine element in hanging node connections.
  Sequence<Element_connection&>& element_connections() {return elem_cons;}
  int add_boundary_condition(Flow_bc*, Mesh_bc*) override;
  void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n) override;
  void disconnect_boundary(int bc_sn) override;

  void add_tree(std::vector<Flow_bc*> extremal_bcs) override;
  void set_surface(Surface_geom* geometry, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start = Eigen::VectorXd::Zero(3)) override;
  void update(std::function<bool(Element&)> refine_criterion = always, std::function<bool(Element&)> unrefine_criterion = never) override;

  //! \returns a view of all Bounday_condition objects owned by this mesh
  Vector_view<Boundary_condition&, Boundary_condition> boundary_conditions() {return bound_conds;}
  //! get a boundary condition owned by this mesh by its serial number
  Boundary_condition& boundary_condition(int bc_sn) {return bound_conds[bc_sn];}
  //! \returns a view of all connections between an element and a boundary condition
  Sequence<Boundary_connection&>& boundary_connections() {return bound_cons;}
  //! \returns a view of all Refined_face objects owned by this mesh (there will be one for every hanging node connection)
  inline Sequence<Refined_face&>& refined_faces() {return ref_face_v;}
  //! \returns a view of all Hanging_vertex_matcher objects owned by this mesh (there will be one for every hanging node connection)
  inline Sequence<Hanging_vertex_matcher&>& hanging_vertex_matchers() {return matcher_v;}
  Connection_validity valid() override;
  //! convenience typedef for the Vector_view used to access Vertex objects
  typedef Vector_view<Vertex&, Vertex::Non_transferable_ptr, &ptr_convert<Vertex&, Vertex::Non_transferable_ptr>> vertex_view;
  //! \returns a view of all Vertex objects used by elements in this mesh. Each vertex will appear exactly once, even if it is shared by multiple elements.
  vertex_view vertices();
  void extrude(bool collapse = false, double offset = 0) override; // note: test for this is in `test_Solver.cpp` so that the result can be visualized
  void connect_rest(int bc_sn) override;
  std::vector<elem_handle> elem_handles() override;
  //! \returns a view of all Element_connection between extruded elements and the elemens they were extruded from
  inline Index<Element_connection&> extruded_connections() {return {deformed().element_connections(), extrude_cons};}
};

}
#endif
