#ifndef HEXED_ELEMENT_CONTAINER_HPP_
#define HEXED_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
#include "math.hpp"
#include "Vector_view.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"

namespace hexed
{

/*!
 * An abstract class representing a collection of `Element`s which supports
 * basic construction and access facilities.
 */
class Element_container
{
  public:
  /*!
   * Construct an element, add it to the container, and return a permanent, arbitrary
   * serial number which is unique among elements of the same refinement level.
   */
  virtual int emplace(int ref_level, std::vector<int> position) = 0;
  //! access an element by refinement level and serial number
  virtual Element& at(int ref_level, int serial_n) = 0;
  virtual Sequence<Element&>& element_view() = 0;
  //! return the currently valid set of `ref_level`, `serial_n` combinations
  virtual std::vector<std::array<int, 2>> elem_handles() = 0;
};

/*!
 * An implementation of `Element_container` which holds elements of a definite type
 * (namely, `Element` or `Deformed_element`). Supports more advanced access methods.
 */
template <typename element_t>
class Complete_element_container : public Element_container
{
  typedef std::unique_ptr<element_t> ptr_t;
  static element_t& convert(ptr_t& ptr) {return *ptr;}
  Storage_params params;
  double spacing;
  int next_sn;
  std::vector<ptr_t> vec;
  std::map<std::array<int, 2>, element_t&> map;
  Vector_view<Element&, ptr_t, ptr_convert<Element&, ptr_t>> view;

  public:
  //! construct with the same data as `Mesh::Mesh`
  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}, next_sn{0}, view{vec}
  {}

  int emplace(int ref_level, std::vector<int> position) override
  {
    vec.emplace_back(new element_t {params, position, spacing, ref_level});
    std::array<int, 2> key = {ref_level, next_sn++};
    map.insert(std::pair<std::array<int, 2>, element_t&>(key, *vec.back()));
    return key[1];
  }

  element_t& at(int ref_level, int serial_n) override
  {
    return map.at({ref_level, serial_n});
  }

  typedef Vector_view<element_t&, ptr_t, &convert> view_t; //!< convenience typedef
  //! Provides a `Vector_view` which can be used to efficiently iterate through the elements, in order of creation (oldest first)
  view_t elements()
  {
    return vec;
  }

  Sequence<Element&>& element_view() {return view;} //!< same as `elements()` except views elements as type `Element&`

  std::vector<std::array<int, 2>> elem_handles() override
  {
    std::vector<std::array<int, 2>> handles;
    for (auto pair : map) handles.push_back(pair.first);
    return handles;
  }
};

}
#endif
