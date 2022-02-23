#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
#include "math.hpp"
#include "Element.hpp"
#include "Deformed_element.hpp"

namespace cartdg
{

class Element_container
{
  public:
  class Mesh_element
  {
    public:
    std::array<int, 6> connectedness;
    virtual Element& get() = 0;
  };
  virtual int emplace(int ref_level, std::vector<int> position) = 0;
  virtual Mesh_element& at(int ref_level, int serial_n) = 0;
};

template <typename element_t>
class Complete_element_container : public Element_container
{
  public:
  class Compl_mesh_elem : public Mesh_element
  {
    public:
    element_t element;
    Compl_mesh_elem(Storage_params par, std::vector<int> pos, double sz) : element{par, pos, sz} {connectedness.fill(0);}
    virtual Element& get() {return element;}
  };

  Complete_element_container(Storage_params storage_params, double root_spacing)
  : params{storage_params}, spacing{root_spacing}, next_sn{0}
  {}

  virtual int emplace(int ref_level, std::vector<int> position)
  {
    double level_spacing = spacing/custom_math::pow(2, ref_level);
    vec.emplace_back(new Compl_mesh_elem {params, position, level_spacing});
    std::array<int, 2> key = {ref_level, next_sn++};
    map.insert(std::pair<std::array<int, 2>, Compl_mesh_elem&>(key, *vec.back()));
    return key[1];
  }

  virtual Compl_mesh_elem& at(int ref_level, int serial_n)
  {
    return map.at({ref_level, serial_n});
  }

  private:
  Storage_params params;
  double spacing;
  int next_sn;
  std::vector<std::unique_ptr<Compl_mesh_elem>> vec;
  std::map<std::array<int, 2>, Compl_mesh_elem&> map;
};

}
#endif
