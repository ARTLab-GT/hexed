#ifndef CARTDG_ACCESSIBLE_MESH_HPP_
#define CARTDG_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"
#include "connection.hpp"

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
  Complete_element_container<Element>          car_elems;
  Complete_element_container<Deformed_element> def_elems;
  Element_container& container(bool is_deformed);
  std::vector<Cartesian_element_connection> car_cons;
  std::vector< Deformed_element_connection> def_cons;

  public:
  class Element_sequence : public Sequence<Element&>
  {
    Accessible_mesh& am;
    public:
    Element_sequence(Accessible_mesh&);
    int size();
    Element& operator[](int);
  };

  Accessible_mesh(Storage_params, double root_size);
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position);
  // Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
  // The following 3 functions return objects which support efficient read/write access to the elements
  // with `size()` and `operator[]` member functions. Elements are not guaranteed to be in any particular order.
  inline auto cartesian_elements() {return car_elems.elements();} // access only cartesian elements
  inline auto  deformed_elements() {return def_elems.elements();} // access only deformed
  Element_sequence elements(); // access all elements regardless of deformedness
  virtual void connect_cartesian(int ref_level, int i_dim, std::array<int, 2> serial_n,
                                 std::array<bool, 2> is_deformed = {false, false});
  // Provides read access to all connections between Cartesian elements in unspecified order
  Vector_view<Cartesian_face_connection&, Cartesian_element_connection> cartesian_connections() {return car_cons;}
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, std::array<int, 2> i_dim, std::array<bool, 2> face_sign);
  // Provides read access to all connections between deformed elements in unspecified order
  Vector_view<Deformed_face_connection&, Deformed_element_connection> deformed_connections() {return def_cons;}
};

}
#endif
