#ifndef CARTDG_ACCESSIBLE_MESH_HPP_
#define CARTDG_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"

namespace cartdg
{

// Provides access to all of the elements, connections, and other numerical data of a specific type
// (i.e. Cartesian or deformed) without addition/removal
template <typename element_t>
class View_by_type
{
  public:
  typedef typename Complete_element_container<element_t>::view_t element_view;
  typedef Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> connection_view;
  virtual element_view elements() = 0;
  virtual connection_view connections() = 0;
};

/*
 * A mesh that supports access to the actual elements with the numerical data they contain. This level of
 * access is required by the numerical scheme but should be hidden from the library user, who should not be
 * concerned with numerical details.
 */
class Accessible_mesh : public Mesh
{
  // Stores numerical data of a particular type and provides free access. Private because it allows
  // open access rather than a neat interface.
  template <typename element_t>
  class Mesh_by_type : public View_by_type<element_t>
  {
    public:
    Complete_element_container<element_t> elems;
    std::vector<Element_face_connection<element_t>> cons;
    Mesh_by_type(Storage_params params, double root_spacing) : elems{params, root_spacing} {}
    virtual typename View_by_type<element_t>::element_view elements() {return elems.elements();}
    virtual typename View_by_type<element_t>::connection_view connections() {return cons;}
  };

  Storage_params params;
  double root_sz;
  Element_container& container(bool is_deformed);
  Mesh_by_type<         Element> car;
  Mesh_by_type<Deformed_element> def;

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
  inline View_by_type<         Element>& cartesian() {return car;}
  inline View_by_type<Deformed_element>&  deformed() {return def;}
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position);
  // Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
  // access all elements regardless of deformedness
  Element_sequence elements();
  virtual void connect_cartesian(int ref_level, std::array<int, 2> serial_n, Con_dir<Element> dir,
                                 std::array<bool, 2> is_deformed = {false, false});
  virtual void connect_deformed(int ref_level, std::array<int, 2> serial_n, Con_dir<Deformed_element> direction);
  #if 0
  virtual void connect_hanging_cartesian(int coarse_ref_level, int coarse_serial, std::vector<int> fine_serial,
                                         Con_dir<Element>, bool coarse_face_positive);
  #endif
};

}
#endif
