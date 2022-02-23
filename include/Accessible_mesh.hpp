#ifndef CARTDG_ACCESSIBLE_MESH_HPP_
#define CARTDG_ACCESSIBLE_MESH_HPP_

#include "Mesh.hpp"
#include "Storage_params.hpp"
#include "Element.hpp"

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

  struct element_handle {int ref_level, int serial_n};
  template<class element_t>
  class element_container {
    public:
    int next_sn = 0;
    std::vector<std::unique_ptr<element_t>> elements;
    std::map<element_handle, element_t*> directory;
    int push(Storage_params sp, std::vector<int> pos, double sz)
    {
      elements.emplace_back(new element_t {sp, pos, sz});
      return next_sn++;
    }
  }
  element_container<Element> cartesian;
  element_container<Deformed_element> deformed;

  public:
  Accessible_mesh(Storage_params, double root_size);
  virtual int add_element(int ref_level, bool is_deformed, std::vector<int> position);

  // Access an element. If the parameters to not describe an existing element, throw an exception.
  Element& element(int ref_level, bool is_deformed, int serial_n);
};

}
#endif
