#ifndef CARTDG_MESH_BY_TYPE_HPP_
#define CARTDG_MESH_BY_TYPE_HPP_

#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"

namespace cartdg
{

/*
 * Provides access to all of the elements, connections, and other numerical data of a specific type
 * (i.e. Cartesian or deformed) without addition/removal. This is a suitable interface for the
 * numerical scheme to interact with.
 */
template <typename element_t>
class View_by_type
{
  public:
  virtual ~View_by_type() = default;
  virtual Sequence<element_t&>& elements() = 0;
  virtual Sequence<Face_connection<element_t>&>& face_connections() = 0;
  virtual Sequence<Element_connection&>& element_connections() = 0;
  virtual Sequence<Refined_face&>& refined_faces() = 0;
  virtual Sequence<Boundary_connection&>& boundary_connections() = 0;
};

/*
 * Stores numerical data of a particular type and provides free access. This is really a helper class
 * for `Accessible_mesh` which grew to the point that it deserved its own file.
 */
template <typename element_t>
class Mesh_by_type : public View_by_type<element_t>
{
  public:
  // where the data is kept
  Complete_element_container<element_t> elems;
  std::vector<Element_face_connection<element_t>> cons;
  std::vector<Refined_connection<element_t>> ref_face_cons;
  std::vector<Typed_bound_connection<element_t>> bound_cons;

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
      return parent.ref_face_cons[i_refined/4].connection(i_refined%4);
    }
  };
  Connection_view<Face_connection<element_t>&> elem_face_con_v;
  // this is useful to allow optional concatenation by providing an empty vector to concatenate
  static std::vector<Element_face_connection<element_t>> empty_con_vec;
  static Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> empty_con_view;
  Concatenation<Face_connection<element_t>&> face_con_v;
  Connection_view<Element_connection&> elem_con_v;
  static Refined_face& ref_face(Refined_connection<element_t>& ref_con) {return ref_con.refined_face;}
  Vector_view<Refined_face&, Refined_connection<element_t>, &ref_face> ref_v;
  Vector_view<Boundary_connection&, Typed_bound_connection<element_t>> bound_con_v;

  public:
  Vector_view<Face_connection<Deformed_element>&, Typed_bound_connection<element_t>> bound_face_con_view;

  // `extra_con_v` let's us tack some extra connections on to the connection view if we want to.
  // In particular, `Accessible_mesh` uses this to put all the boundary connections into the deformed connections
  Mesh_by_type(Storage_params params, double root_spacing,
               Sequence<Face_connection<element_t>&>& extra_con_v = empty_con_view)
  : elems{params, root_spacing},
    elem_v{elems.elements()},
    elem_face_con_v{*this},
    face_con_v{elem_face_con_v, extra_con_v},
    elem_con_v{*this},
    ref_v{ref_face_cons},
    bound_con_v{bound_cons},
    bound_face_con_view{bound_cons}
  {}

  // `View_by_type` interface implementation
  virtual Sequence<element_t&>& elements() {return elem_v;}
  virtual Sequence<Face_connection<element_t>&>& face_connections() {return face_con_v;}
  virtual Sequence<Element_connection&>& element_connections() {return elem_con_v;}
  virtual Sequence<Refined_face&>& refined_faces() {return ref_v;}
  virtual Sequence<Boundary_connection&>& boundary_connections() {return bound_con_v;}
};

template <typename element_t>
std::vector<Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_vec {};

template <typename element_t>
Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_view {Mesh_by_type<element_t>::empty_con_vec};

}
#endif
