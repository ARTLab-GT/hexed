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
  //!< deletes all elements where `predicate` evaluates to true and returns the number of elements deleted
  virtual int purge(std::function<bool(Element&)> predicate = [](Element& elem){return elem.record != 0;}) = 0;
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
  std::vector<int> serial_ns;
  Vector_view<Element&, ptr_t, ptr_convert<Element&, ptr_t>> view;

  public:
  //! construct with the same data as `Mesh::Mesh`
  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}, next_sn{0}, view{vec}
  {}

  int emplace(int ref_level, std::vector<int> position) override
  {
    vec.emplace_back(new element_t {params, position, spacing, ref_level});
    serial_ns.push_back(next_sn++);
    return serial_ns.back();
  }

  element_t& at(int ref_level, int serial_n) override
  {
    auto range = std::equal_range(serial_ns.begin(), serial_ns.end(), serial_n);
    HEXED_ASSERT(range.second - range.first == 1, "could not find unique element with specified serial number");
    element_t& elem = *vec[range.first - serial_ns.begin()];
    HEXED_ASSERT(elem.refinement_level() == ref_level, "could not find element with specified refinement level");
    return elem;
  }

  typedef Vector_view<element_t&, ptr_t, &convert> view_t; //!< convenience typedef
  //! Provides a `Vector_view` which can be used to efficiently iterate through the elements, in order of creation (oldest first)
  view_t elements()
  {
    return vec;
  }

  Sequence<Element&>& element_view() {return view;} //!< same as `elements()` except views elements as type `Element&`

  int purge(std::function<bool(Element&)> predicate = [](Element& elem){return elem.record != 0;}) override
  {
    int old_size = vec.size();
    std::vector<std::unique_ptr<element_t>> new_vec;
    std::vector<int> new_sns;
    for (int i_elem = 0; i_elem < old_size; ++i_elem) {
      if (predicate(*vec[i_elem])) {
        if (vec[i_elem]->tree) vec[i_elem]->tree->elem = nullptr;
      } else {
        new_vec.emplace_back(vec[i_elem].release());
        new_sns.push_back(serial_ns[i_elem]);
      }
    }
    vec = std::move(new_vec);
    serial_ns = std::move(new_sns);
    return old_size - vec.size();
  }

  std::vector<std::array<int, 2>> elem_handles() override
  {
    std::vector<std::array<int, 2>> handles;
    for (unsigned i_elem = 0; i_elem < vec.size(); ++i_elem) {
      handles.push_back({vec[i_elem]->refinement_level(), serial_ns[i_elem]});
    }
    return handles;
  }

  //! writes each element's serial number to it's `record` \todo elements really just need to know their serial number
  void write_sns()
  {
    #pragma omp parallel for
    for (unsigned i_elem = 0; i_elem < vec.size(); ++i_elem) {
      vec[i_elem]->record = serial_ns[i_elem];
    }
  }
};

}
#endif
