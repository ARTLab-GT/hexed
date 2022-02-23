#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
#include "math.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

class Element_vector
{
  public:
  virtual int size() = 0;
  virtual Element& operator[](int index) = 0;
};

class Element_container : public Element_vector
{
  virtual int emplace(int ref_level, std::vector<int> position) = 0;
  virtual Element& at(int ref_level, int serial_n) = 0;
  virtual std::array<int, 6> connectedness(int ref_level, int serial_n) = 0;
};

template <typename element_t>
class Complete_element_container : public Element_container
{
  Storage_params params;
  double spacing;
  int next_sn;
  std::vector<std::unique_ptr<element_t>> vec;
  struct Element_data {element_t& element; std::array<int, 6> connectedness;};
  std::map<std::array<int, 2>, Element_data> map;

  public:
  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}, next_sn{0}
  {}
  virtual int size()
  {
    return 0;
  }
  virtual element_t& operator[](int index)
  {
    return *vec.back();
  }
  virtual int emplace(int ref_level, std::vector<int> position)
  {
    double level_spacing = spacing/custom_math::pow(2, ref_level);
    vec.emplace_back(new element_t {params, position, level_spacing});
    std::array<int, 2> key = {ref_level, next_sn++};
    Element_data data {*vec.back(), {}};
    map.insert(std::pair(key, data));
    return key[1];
  }
  virtual element_t& at(int ref_level, int serial_n)
  {
    return map.at({ref_level, serial_n}).element;
  }
  virtual std::array<int, 6> connectedness(int ref_level, int serial_n)
  {
    return {0, 0, 0, 0, 0, 0};
  }
};

}
#endif
