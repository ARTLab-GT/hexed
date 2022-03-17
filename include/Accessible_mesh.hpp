#ifndef CARTDG_ACCESSIBLE_MESH_HPP_
#define CARTDG_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"

namespace cartdg
{

/*
 * Provides access to all of the elements, connections, and other numerical data of a specific type
 * (i.e. Cartesian or deformed) without addition/removal
 */
template <typename element_t>
class View_by_type
{
  public:
  virtual Sequence<element_t&>& elements() = 0;
  virtual Sequence<Face_connection<element_t>&>& face_connections() = 0;
  virtual Sequence<Element_connection&>& element_connections() = 0;
  virtual Sequence<Refined_face&>& refined_faces() = 0;
};

// Stores numerical data of a particular type and provides free access.
template <typename element_t>
class Mesh_by_type : public View_by_type<element_t>
{
  public:
  // where the data is kept
  Complete_element_container<element_t> elems;
  std::vector<Element_face_connection<element_t>> cons;
  std::vector<Refined_connection<element_t>> ref_face_cons;
  private:
  // template spaghetti to get `Vector_view`s of the data with the right type
  typename Complete_element_container<element_t>::view_t elem_v;
  template <typename view_t>
  class Connection_view : public Sequence<view_t>
  {
    Mesh_by_type& parent;
    public:
    Connection_view(Mesh_by_type& mbt) : parent{mbt} {}
    virtual int size() {return parent.cons.size() + 4*parent.ref_face_cons.size();}
    virtual view_t operator[](int index)
    {
      int i_refined = index - parent.cons.size();
      if (i_refined < 0) return parent.cons[index];
      else return parent.ref_face_cons[i_refined/4].connection(i_refined%4);
    }
  };
  Connection_view<Face_connection<element_t>&> face_con_v;
  Connection_view<Element_connection&> elem_con_v;
  static Refined_face& ref_face(Refined_connection<element_t>& ref_con) {return ref_con.refined_face;}
  Vector_view<Refined_face&, Refined_connection<element_t>, &ref_face> ref_v;
  public:
  Mesh_by_type(Storage_params params, double root_spacing)
  : elems{params, root_spacing}, elem_v{elems.elements()}, face_con_v{*this},
    elem_con_v{*this}, ref_v{ref_face_cons}
  {}
  // interface implementation
  virtual Sequence<element_t&>& elements() {return elem_v;}
  virtual Sequence<Face_connection<element_t>&>& face_connections() {return face_con_v;}
  virtual Sequence<Element_connection&>& element_connections() {return elem_con_v;}
  virtual Sequence<Refined_face&>& refined_faces() {return ref_v;}
};

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
};

}
#endif
