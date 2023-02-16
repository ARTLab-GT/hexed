#ifndef HEXED_ACCESSIBLE_MESH_HPP_
#define HEXED_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Mesh_by_type.hpp"

namespace hexed
{

/*
 * A mesh that supports access to the actual elements with the numerical data they contain. This level of
 * access is required by the numerical scheme but should be hidden from the library user, who should not be
 * concerned with numerical details.
 */
class Accessible_mesh : public Mesh
{
  Storage_params params;
  int n_vert;
  double root_sz;
  Element_container& container(bool is_deformed);
  Mesh_by_type<         Element> car;
  Mesh_by_type<Deformed_element> def;
  // create a `Vector_view` that can look at `def.elements()` as `Element&`s.
  Vector_view<Element&, Deformed_element&, &trivial_convert<Element&, Deformed_element&>, Sequence> def_as_car;
  Concatenation<Element&> elems;
  Concatenation<Element_connection&> elem_cons;
  std::vector<Boundary_condition> bound_conds;
  Concatenation<Face_connection<Deformed_element>&> bound_face_cons;
  Concatenation<Boundary_connection&> bound_cons;
  Concatenation<Face_connection<Deformed_element>&> def_face_cons;
  Concatenation<Refined_face&> ref_face_v;
  Concatenation<Hanging_vertex_matcher&> matcher_v;
  std::vector<Vertex::Non_transferable_ptr> vert_ptrs;
  std::vector<int> extrude_cons;

  public:
  Accessible_mesh(Storage_params, double root_size);
  virtual inline double root_size() {return root_sz;}
  inline View_by_type<         Element>& cartesian() {return car;}
  inline View_by_type<Deformed_element>&  deformed() {return def;}
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position);
  // Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
  // access all elements regardless of deformedness
  Sequence<Element&>& elements() {return elems;}
  virtual void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                                 std::array<bool, 2> is_deformed = {false, false});
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction);
  virtual void connect_hanging(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Deformed_element>,
                               bool coarse_deformed = false, std::vector<bool> fine_deformed = {false, false, false, false},
                               std::array<bool, 2> stretch = {false, false});
  Sequence<Element_connection&>& element_connections() {return elem_cons;}
  virtual int add_boundary_condition(Flow_bc*, Mesh_bc*);
  virtual void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n);
  virtual void disconnect_boundary(int bc_sn);
  Vector_view<Boundary_condition&, Boundary_condition> boundary_conditions() {return bound_conds;}
  Boundary_condition& boundary_condition(int bc_sn) {return bound_conds[bc_sn];}
  Sequence<Boundary_connection&>& boundary_connections() {return bound_cons;}
  inline Sequence<Refined_face&>& refined_faces() {return ref_face_v;}
  inline Sequence<Hanging_vertex_matcher&>& hanging_vertex_matchers() {return matcher_v;}
  virtual Connection_validity valid();
  typedef Vector_view<Vertex&, Vertex::Non_transferable_ptr, &ptr_convert<Vertex&, Vertex::Non_transferable_ptr>> vertex_view;
  vertex_view vertices();
  virtual void extrude(bool collapse = false, double offset = 0); // note: test for this is in `test_Solver.cpp` so that the result can be visualized
  virtual void connect_rest(int bc_sn);
  virtual std::vector<elem_handle> elem_handles();
  inline Index<Element_connection&> extruded_connections() {return {deformed().element_connections(), extrude_cons};}
};

}
#endif
