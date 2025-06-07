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
  //!\}

  /*! \name views
   * template spaghetti to get `Sequence`s of the data with the right type
   */
  //!\{
  typename Complete_element_container<element_t>::view_t elem_v; //!< a view of the elements that does not allow addition or removal
  Vector_view<Kernel_element&, element_t&, &trivial_convert<Kernel_element&, element_t&>, Sequence> kernel_elems;
  //!\}

  Mesh_by_type(Storage_params params, double root_spacing)
  : n_faces{2*params.n_dim}
  , par{params}
  , elems{params, root_spacing}
  , elem_v{elems.elements()}
  , kernel_elems{elem_v}
  {}

  // `View_by_type` interface implementation
  Sequence<element_t&>& elements() override {return elem_v;}
  Sequence<Kernel_element&>& kernel_elements() override {return kernel_elems;}

};

}
#endif
