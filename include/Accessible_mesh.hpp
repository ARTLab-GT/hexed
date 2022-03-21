#ifndef CARTDG_ACCESSIBLE_MESH_HPP_
#define CARTDG_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Mesh_by_type.hpp"

namespace cartdg
{

/*
 * A mesh that supports access to the actual elements with the numerical data they contain. This level of
 * access is required by the numerical scheme but should be hidden from the library user, who should not be
 * concerned with numerical details.
 */
class Accessible_mesh : public Mesh
{
  Storage_params params;
  double root_sz;
  Element_container& container(bool is_deformed);
  Mesh_by_type<         Element> car;
  Mesh_by_type<Deformed_element> def;
  // create a `Vector_view` that can look at `def.elements()` as `Element&`s.
  Vector_view<Element&, Deformed_element&, &trivial_convert<Element&, Deformed_element&>, Sequence> def_as_car;
  Concatenation<Element&> elems;
  Concatenation<Element_connection&> elem_cons;
  std::vector<std::unique_ptr<Boundary_condition>> bound_conds;
  static Boundary_condition& convert_bc (std::unique_ptr<Boundary_condition>& ptr) {return *ptr;}
  Vector_view<Boundary_condition&, std::unique_ptr<Boundary_condition>, &convert_bc> bc_v;
  Concatenation<Face_connection<Deformed_element>&> bound_face_cons;
  Concatenation<Face_connection<Deformed_element>&> def_face_cons;

  public:
  Accessible_mesh(Storage_params, double root_size);
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
  virtual void connect_hanging_cartesian(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial, Con_dir<Element>,
                                         bool coarse_face_positive, bool coarse_deformed=false, bool fine_deformed=false);
  Sequence<Element_connection&>& element_connections() {return elem_cons;}
  virtual int add_boundary_condition(Boundary_condition* bc);
  virtual void connect_boundary(int ref_level, bool is_deformed, int element_serial_n, int i_dim, int face_sign, int bc_serial_n);
  Sequence<Boundary_condition&>& boundary_conditions() {return bc_v;}
  virtual Connection_validity valid();
};

}
#endif
