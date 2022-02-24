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
 * A mesh that supports access to the actual elements with the numerical data they contain. This level of
 * access is required by the numerical scheme but should be hidden from the library user, which should not be
 * concerned with numerical details.
 */
class Accessible_mesh : public Mesh
{
  Storage_params params;
  double root_sz;
  Complete_element_container<Element>          car_elements;
  Complete_element_container<Deformed_element> def_elements;
  Element_container& container(bool is_deformed);

  public:
  Accessible_mesh(Storage_params, double root_size);
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position);
  // Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
};

}
#endif
