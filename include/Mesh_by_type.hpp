#ifndef HEXED_MESH_BY_TYPE_HPP_
#define HEXED_MESH_BY_TYPE_HPP_

#include "Storage_params.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"
#include "Element_container.hpp"
#include "Tree.hpp"
#include "utils.hpp"

namespace hexed
{

/*!
 * Provides access to all of the elements, connections, and other numerical data of a specific type
 * (i.e. Cartesian or deformed) without addition/removal. This is a suitable interface for the
 * numerical scheme to interact with.
 */
template <typename element_t>
class View_by_type
{
  public:
  virtual ~View_by_type() = default;
  //! \cond
  virtual Sequence<element_t&>& elements() = 0;
  virtual Sequence<Kernel_element&>& kernel_elements() = 0;
  virtual Sequence<Face_connection<element_t>&>& face_connections() = 0;
  virtual Sequence<Kernel_connection&>& kernel_connections() = 0;
  virtual Sequence<Element_connection&>& element_connections() = 0;
  //! \endcond
};

/*!
 * Stores numerical data of a particular type and provides free access.
 * This is really a helper class for `Accessible_mesh`
 * which grew to the point that it deserved its own file.
 */
template <typename element_t>
class Mesh_by_type : public View_by_type<element_t>
{
  const int n_faces;
  Storage_params par;

  public:
  /*! \name containers
   * where the actual data is kept
   */
  //!\{
  Complete_element_container<element_t> elems;
  std::vector<std::unique_ptr<Element_face_connection<element_t>>> cons;
  //!\}

  /*! \name views
   * template spaghetti to get `Sequence`s of the data with the right type
   */
  //!\{
  typename Complete_element_container<element_t>::view_t elem_v; //!< a view of the elements that does not allow addition or removal
  Vector_view<Kernel_element&, element_t&, &trivial_convert<Kernel_element&, element_t&>, Sequence> kernel_elems;

  //! Sequence of some type of connection object which cycles through first the conformal connections and then the hanging-node connections.
  template <typename view_t>
  class Connection_view : public Sequence<view_t>
  {
    Mesh_by_type& parent;
    public:
    Connection_view(Mesh_by_type& mbt) : parent{mbt} {}
    virtual int size()
    {
      int sz = parent.cons.size();
      return sz;
    }
    virtual view_t operator[](int index)
    {
      int i_start = index - parent.cons.size();
      if (i_start < 0) return *parent.cons[index];
      throw std::runtime_error("`Connection_view` indexed out of bounds.");
    }
  };
  Connection_view<Face_connection<element_t>&> elem_face_con_v;
  // this is useful to allow optional concatenation by providing an empty vector to concatenate
  static std::vector<Element_face_connection<element_t>> empty_con_vec;
  static Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> empty_con_view;
  Connection_view<Element_connection&> elem_con_v;
  Vector_view<Kernel_connection&, Face_connection<element_t>&, &trivial_convert<Kernel_connection&, Face_connection<element_t>&>, Sequence> kernel_cons;
  //!\}

  Mesh_by_type(Storage_params params, double root_spacing)
  : n_faces{2*params.n_dim}
  , par{params}
  , elems{params, root_spacing}
  , elem_v{elems.elements()}
  , kernel_elems{elem_v}
  , elem_face_con_v{*this}
  , elem_con_v{*this}
  , kernel_cons{elem_face_con_v}
  {}

  // `View_by_type` interface implementation
  Sequence<element_t&>& elements() override {return elem_v;}
  Sequence<Kernel_element&>& kernel_elements() override {return kernel_elems;}
  Sequence<Face_connection<element_t>&>& face_connections() override {return elem_face_con_v;}
  Sequence<Kernel_connection&>& kernel_connections() override {return kernel_cons;}
  Sequence<Element_connection&>& element_connections() override {return elem_con_v;}

  //! write the number of connections for each face to `Element::face_record`. Assumes initialized to 0
  void record_connections()
  {
    // ordinary element connections
    for (unsigned i_con = 0; i_con < cons.size(); ++i_con) {
      for (int i_side = 0; i_side < 2; ++i_side) {
        ++cons[i_con]->element(i_side).face_record[cons[i_con]->direction().i_face(i_side)];
      }
    }
  }

  //! delete all connections (of all kinds) where `predicate` is true for at least one of the elements involved
  //! or connections for boundary faces that have since been covered up
  void purge_connections(std::function<bool(Element&)> predicate = [](Element& elem){return elem.record != 0;})
  {
    erase_if(cons, [predicate](std::unique_ptr<Element_face_connection<element_t>>& con) {
      if (!con) return true;
      return predicate(con->element(0)) || predicate(con->element(1)) || !con->neighbor_connection().alive();
    });
  }
};

template <typename element_t>
std::vector<Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_vec {};

template <typename element_t>
Vector_view<Face_connection<element_t>&, Element_face_connection<element_t>> Mesh_by_type<element_t>::empty_con_view {Mesh_by_type<element_t>::empty_con_vec};

}
#endif
