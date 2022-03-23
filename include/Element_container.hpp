#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
#include "math.hpp"
#include "Vector_view.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

/*
 * An abstract class representing a collection of `Element`s which supports
 * basic construction and access facilities.
 */
class Element_container
{
  public:
  /*
   * Construct an element, add it to the container, and return a permanent, arbitrary
   * serial number which is unique among elements of the same refinement level.
   */
  virtual int emplace(int ref_level, std::vector<int> position) = 0;
  // access an element by refinement level and serial number
  virtual Element& at(int ref_level, int serial_n) = 0;
};

/*
 * An implementation of `Element_container` which holds elements of a definite type (namely,
 * `Element` or `Deformed_element`). Supports more advanced access methods.
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

  public:
  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}, next_sn{0}
  {}

  virtual int emplace(int ref_level, std::vector<int> position)
  {
    double level_spacing = spacing/custom_math::pow(2, ref_level);
    vec.emplace_back(new element_t {params, position, level_spacing});
    std::array<int, 2> key = {ref_level, next_sn++};
    map.insert(std::pair<std::array<int, 2>, element_t&>(key, *vec.back()));
    return key[1];
  }

  virtual element_t& at(int ref_level, int serial_n)
  {
    return map.at({ref_level, serial_n});
  }

  // Provides a `Vector_view` which can be used to efficiently iterate through the elements,
  // in no particular order.
  typedef Vector_view<element_t&, ptr_t, &convert> view_t;
  view_t elements()
  {
    return vec;
  }
};

}
#endif
