#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
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
  virtual void erase(int ref_level, int serial_n) = 0;
  virtual Element& at(int ref_level, int serial_n) = 0;
  virtual std::array<int, 6> connectedness(int ref_level, int serial_n) = 0;
};

template <typename element_t>
class Complete_element_container : public Element_container
{
  Storage_params params;
  double spacing;
  std::vector<std::unique_ptr<element_t>> vec;
  struct Element_data {element_t& element; std::array<int, 6> connectedness;};
  std::map<std::array<int, 2>, Element_data> map;

  public:
  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}
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
    return 0;
  }
  virtual void erase(int ref_level, int serial_n)
  {
  }
  virtual element_t& at(int ref_level, int serial_n)
  {
    return *vec.back();
  }
  virtual std::array<int, 6> connectedness(int ref_level, int serial_n)
  {
    return {0, 0, 0, 0, 0, 0};
  }
};

}
#endif
