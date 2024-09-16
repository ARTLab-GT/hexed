#ifndef HEXED_ACCESSIBLE_MESH_HPP_
#define HEXED_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Mesh_by_type.hpp"
#include "Tree.hpp"
#include "Kernel_mesh.hpp"
#include "Gauss_lobatto.hpp"

namespace hexed {

/*! \brief A mesh that supports access to the actual elements with the numerical data they contain.
 * \details This level of access is required by the numerical scheme but should be hidden from the library user,
 * who should not be concerned with numerical details.
 */
class Accessible_mesh : public Mesh {
  Storage_params params;
  int n_vert;
  double root_sz;
  Mesh_by_type<         Element> car;
  Mesh_by_type<Deformed_element> def;
  // create a `Vector_view` that can look at `def.elements()` as `Element&`s.
  Vector_view<Element&, Deformed_element&, &trivial_convert<Element&, Deformed_element&>, Sequence> def_as_car;
  Concatenation<Element&> elems;
  Vector_view<Kernel_element&, Element&, &trivial_convert<Kernel_element&, Element&>, Sequence> kernel_elems;
  Concatenation<Element_connection&> elem_cons;
  std::vector<std::unique_ptr<Flow_bc>> bound_conds;
  Concatenation<Face_connection<Deformed_element>&> bound_face_cons;
  Concatenation<Boundary_connection&> bound_cons;
  Concatenation<Face_connection<Deformed_element>&> def_face_cons;
  Concatenation<Refined_face&> ref_face_v;
  Concatenation<Hanging_vertex_matcher&> matcher_v;
  std::vector<Vertex::Non_transferable_ptr> vert_ptrs;
  int surf_bc_sn;
  std::unique_ptr<Surface_geom> surf_geom;
  std::vector<Element_face_connection<Deformed_element>*> extrude_cons;
  std::unique_ptr<Tree> tree; // could be null! don't forget to check
  std::vector<int> tree_bcs;
  bool verts_are_reset;
  std::vector<std::vector<Vertex::Non_transferable_ptr>> boundary_verts; // a vector of the vertices that are on each boundary
  std::vector<Vertex::Non_transferable_ptr> smooth_verts; // a vector of the vertices that need to be smoothed in this sweep
  int _mask_levels;
  Gauss_lobatto _basis;
  next::Mesh_blocks _blocks;
  int _n_verts;
  Stopwatch_tree _stopwatch;
  std::vector<Mortal_ptr<next::Vertex>> point_matched_vertices;
  std::vector<std::vector<Mortal_ptr<next::Vertex>>> matched_vertices;
  std::vector<std::vector<Mortal_ptr<next::Edge>>> matched_edges;

  // masked sequences
  template <typename view_t, typename storage_t>
  struct Masked {
    std::vector<storage_t*> ptrs;
    Vector_view<view_t&, storage_t*, ptr_convert<view_t&, storage_t*>> view;
    Slice<view_t&> slice;
    Masked() : view(ptrs), slice(view) {}

    template <typename T, typename U>
    void populate(T& base_seq, U criterion) {
      ptrs.resize(base_seq.size());
      int i_ptr = 0;
      for (int i = 0; i < int(base_seq.size()); ++i) {
        if (criterion(base_seq[i])) ptrs[i_ptr++] = &base_seq[i];
      }
      slice = Slice<view_t&>(view, 0, i_ptr);
    }
  };

  Element_container& container(bool is_deformed);
  int add_element(int ref_level, bool is_deformed, std::vector<int> position, Mat<> origin,
                  int aniso_ref_level = 0, int surface_face = next::Mesh_blocks::no_face);
  Element& add_elem(bool is_deformed, Tree&);
  bool intersects_surface(Tree*);
  bool is_surface(Tree*);
  template<typename element_t> Mesh_by_type<element_t>& mbt(); // gets either `car` or `def`
  template<typename element_t> void connect_new(int start_at); // connects new elements in `mbt<element_t>()`. helper function for `refine`
  void refine_set_status(Tree*); // refines a tree and sets the flood fill status for any children that intersect the surface
  void refine_by_record(bool is_deformed, int start, int end);
  bool needs_refine(Tree*);
  void purge();
  void delete_bad_extrusions();
  void deform();
  void id_smooth_verts();
  void id_boundary_verts();
  // identify which vertices are on which boundaries and write it to `Vertex::record`
  // must be called directly befor `snap_vertices`
  void snap_vertices();
  void create_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin = Mat<>::Zero(3));
  void read_file(std::string file_name);
  void _connect(std::array<Element*, 2>, Con_dir<Element>);
  void _connect(std::array<Deformed_element*, 2>, Con_dir<Deformed_element>);
  void _connect(Element*, std::vector<Element*>, Con_dir<Deformed_element>);
  void _connect(Deformed_element*, std::vector<Deformed_element*>, Con_dir<Deformed_element>,
                std::array<bool, 2> = {false, false});

  template <typename Elem_t>
  void _connect_shapes(Elem_t*, std::vector<Elem_t*>, Con_dir<Deformed_element>, std::array<bool, 2>);

  struct Edge_match {
    Mortal_ptr<next::Edge> edge;
    std::array<Geom_edge::Node, 2> nodes;
  };

  public:
  //! \brief how far must the center of an element be from the geometry relative to the nominal size
  //! \details Defaults to \f$ \frac{\sqrt{n_d}}{2} \f$. You can modify it, but it cannot be less than this.
  double buffer_dist;
  /*!
   * \param params parameters specifying what data is stored in each Element (row size, number of dimensions, etc.)
   * \param root_size defines the \ref root_size of the mesh.
   */
  Accessible_mesh(Storage_params params, double root_size);
  /*! \brief Reads mesh from a file created by `Mesh::write`.
   * \details Acquires ownership of boundary condition pointers.
   * This variant is only for tree meshing.
   * The surface boundary condition and geometry arguments must be specified iff the original mesh had a surface geometry (else exception).
   */
  Accessible_mesh(std::string file_name, std::vector<Flow_bc*> extremal_bcs, Surface_geom* = nullptr, Flow_bc* surface_bc = nullptr);
  /*! \brief Reads mesh from a file created by `Mesh::write`.
   * \details Acquires ownership of boundary condition pointers.
   * This variant is not for tree meshing.
   */
  Accessible_mesh(std::string file_name, std::vector<Flow_bc*>);
  virtual ~Accessible_mesh();
  inline double root_size() override {return root_sz;}
  inline Storage_params storage_params() {return params;}
  //! \returns a View_by_type containing only the Cartesian elements in the mesh
  inline View_by_type<         Element>& cartesian() {return car;}
  //! \returns a View_by_type containing only the deformed elements in the mesh
  inline View_by_type<Deformed_element>&  deformed() {return def;}
  int add_element(int ref_level, bool is_deformed, std::vector<int> position) override;
  //! Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
  //! access all elements, both Cartesian and deformed
  Sequence<Element&>& elements() {return elems;}
  Sequence<Kernel_element&>& kernel_elements() {return kernel_elems;}
  void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                         std::array<bool, 2> is_deformed = {false, false}) override;
  void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction) override;
  void connect_hanging(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Deformed_element>,
                       bool coarse_deformed = false, std::vector<bool> fine_deformed = {false, false, false, false},
                       std::array<bool, 2> stretch = {false, false}) override;
  //! \returns a view of all connections between elements, including one connection for every fine element in hanging node connections.
  Sequence<Element_connection&>& element_connections() {return elem_cons;}
  int add_boundary_condition(Flow_bc*) override;
  void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n) override;
  void disconnect_boundary(int bc_sn) override;
  void cleanup() override;

  void add_tree(std::vector<Flow_bc*> extremal_bcs, Mat<> origin = Mat<>::Zero(3)) override;
  void set_surface(Surface_geom* geometry, Flow_bc* surface_bc, Eigen::VectorXd flood_fill_start = Eigen::VectorXd::Zero(3)) override;
  void relax_and_match(int n_relax = 0, double factor = .9) override;
  void set_unref_locks(std::function<bool(Element&)> lock_if = criteria::never) override;
  bool update(std::function<bool(Element&)> refine_criterion = criteria::always, std::function<bool(Element&)> unrefine_criterion = criteria::never) override;
  void set_all_smooth() override;
  void relax(double factor = 0.9) override;
  inline int surface_bc_sn() override {return surf_bc_sn;}
  inline Surface_geom& surface_geometry() {return *surf_geom;}

  /*! \brief Creates a mask that allows kernel operations to be performed on a subset of the elements.
   * \details Supply a function that returns `true` for elements that should be operated on.
   * The masked mesh is described by two objects:
   * - a `Kernel_mesh` which includes all the elements where the mask is `true`
   *   and all the connections and refined faces where the mask is `true` for at least one of the participating elements
   * - a `Sequence` of `Boundary_connection`s which includes all boundary connections for whose element the mask is `true`
   *
   * The members `elements()`, `element_connections()`, etc. will not be affected.
   * This feature has some quirks.
   * The first mask you create will include all elements, regardless of what mask function you supply.
   * After that, each mask will be a subset of all previous masks, again regardless of the mask function.
   * Calling `reset_masks()` invalidates all previous masks and makes it as if you had not yet created any masks
   * (so now your first one will include all elements, etc.).
   * Modifying the mesh invalidates all existing masks but does __not__ perform a reset.
   * After modifying the mesh, you should call `reset_masks()` before making any new masks.
   * \todo Make this more intuitive and less error-prone.
   */
  class Masked_mesh {
    Masked<Kernel_element, Element> _masked_elems;
    Masked<Kernel_element, Element> _masked_car_elems;
    Masked<Kernel_element, Element> _masked_def_elems;
    Masked<Kernel_connection, Kernel_connection> _masked_car_cons;
    Masked<Kernel_connection, Kernel_connection> _masked_def_cons;
    Masked<Refined_face, Refined_face> _masked_ref_faces;
    Masked<Boundary_connection, Boundary_connection> _masked_bound_cons;
    public:
    Masked_mesh(Accessible_mesh&, const Basis&, std::function<bool(Element&)> = [](Element&){return true;});
    Kernel_mesh kernel_mesh;
    Sequence<Boundary_connection&>& bound_cons;
  };
  //! \brief Resets effective number of masks created to 0 and invalidates existing masks.
  //! \see `Masked_mesh`
  void reset_masks();
  /*! \brief Experimental feature. Ignore for now.
   * \details Has to do with an experimental performance-enhancing feature where anisotropic elements are updated more frequently than isotropic ones.
   * Not ready for production use, although it was the motivation for the `Masked_mesh` feature.
   */
  std::vector<std::unique_ptr<Masked_mesh>> preti_masks(const Basis&);

  //! \returns a view of all Bounday_condition objects owned by this mesh
  Vector_view<Flow_bc&, std::unique_ptr<Flow_bc>, &ptr_convert<Flow_bc&, std::unique_ptr<Flow_bc>>>
  boundary_conditions() {return bound_conds;}
  //! get a boundary condition owned by this mesh by its serial number
  Flow_bc& boundary_condition(int bc_sn) {return *bound_conds[bc_sn];}
  //! \returns a view of all connections between an element and a boundary condition
  Sequence<Boundary_connection&>& boundary_connections() {return bound_cons;}
  //! \returns a view of all Refined_face objects owned by this mesh (there will be one for every hanging node connection)
  inline Sequence<Refined_face&>& refined_faces() {return ref_face_v;}
  //! \returns a view of all Hanging_vertex_matcher objects owned by this mesh (there will be one for every hanging node connection)
  inline Sequence<Hanging_vertex_matcher&>& hanging_vertex_matchers() {return matcher_v;}
  inline int n_elements() override {return elements().size();}
  Connection_validity valid() override;
  //! convenience typedef for the Vector_view used to access Vertex objects
  typedef Vector_view<Vertex&, Vertex::Non_transferable_ptr, &ptr_convert<Vertex&, Vertex::Non_transferable_ptr>> vertex_view;
  //! \returns a view of all Vertex objects used by elements in this mesh. Each vertex will appear exactly once, even if it is shared by multiple elements.
  vertex_view vertices();
  void extrude(bool collapse = false, double offset = 0, bool force = false) override; // note: test for this is in `test_Solver.cpp` so that the result can be visualized
  void connect_rest(int bc_sn) override;
  std::vector<elem_handle> elem_handles() override;
  //! \returns a view of all Element_connection between extruded elements and the elemens they were extruded from
  inline Vector_view<Element_connection&, Element_face_connection<Deformed_element>*,
                     ptr_convert<Element_connection&, Element_face_connection<Deformed_element>*>> extruded_connections() {return {extrude_cons};}
  void write(std::string file_name) override;
  void export_polymesh(std::string dir_name) override;
  void visualize(std::string format, std::string file_name) override;
  inline const Stopwatch_tree& stopwatch_tree() const override {return _stopwatch;}

  protected:
  void reset_verts() override;
  void restore_verts() override;
};

}
#endif
