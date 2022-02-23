#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include "Element.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

class Element_container
{
  protected:
  Storage_params params;
  double spacing;

  public:
  inline Element_container(Storage_params storage_params, double root_spacing) : params{storage_params}, spacing{root_spacing} {}
  Element_container(const Element_container&) = delete;
  virtual ~Element_container() = default;
  virtual int add_element(int ref_level, std::vector<int> position) = 0;
  virtual void erase_element(int ref_level, int serial_n) = 0;
};

class Accessible_element_container : public Element_container
{
  public:
  inline Accessible_element_container(Storage_params storage_params, double root_spacing) : Element_container(storage_params, root_spacing) {}
  virtual Element& element(int ref_level, int serial_n) = 0;
};

template <typename element_t>
class Specific_container : public Accessible_element_container
{
  std::vector<std::unique_ptr<element_t>> elems;

  public:
  inline Specific_container(Storage_params storage_params, double root_spacing) : Accessible_element_container(storage_params, root_spacing) {}
  virtual int add_element(int ref_level, std::vector<int> position) {return 0;}
  virtual void erase_element(int ref_level, int serial_n) {}
  virtual Element& element(int ref_level, int serial_n) {return *elems[0];}
  int size() {return 0;}
  element_t& operator[](int index) {return *elems[0];}
};

}
#endif
