#ifndef CARTDG_ELEMENT_CONTAINER_HPP_
#define CARTDG_ELEMENT_CONTAINER_HPP_

#include <vector>
#include <map>
#include "math.hpp"
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
  // contains an `Element` with some additional data relevant to the mesh
  class Mesh_element
  {
    public:
    std::array<int, 6> connectedness;
    virtual Element& get() = 0;
  };
  /*
   * Construct an element, add it to the container, and return a permanent, arbitrary
   * serial number which is unique among elements of the same refinement level.
   */
  virtual int emplace(int ref_level, std::vector<int> position) = 0;
  // access an element by refinement level and serial number
  virtual Mesh_element& at(int ref_level, int serial_n) = 0;
};

/*
 * An implementation of `Element_container` which holds elements of a definite type (namely,
 * `Element` or `Deformed_element`). Supports more advanced access methods.
 */
template <typename element_t>
class Complete_element_container : public Element_container
{
  public:
  class Compl_mesh_elem : public Mesh_element
  {
    public:
    element_t element;
    Compl_mesh_elem(Storage_params par, std::vector<int> pos, double sz) : element{par, pos, sz} {connectedness.fill(0);}
    operator element_t&() {return element;}
    virtual Element& get() {return element;}
  };
  /*
   * A sequence-style object which provides efficient access to the elements of a
   * `Complete_element_container`, in no particular order.
   */
  class Sequence
  {
    Complete_element_container& c;
    public:
    Sequence(Complete_element_container& container) : c{container} {}
    int size() {return c.vec.size();}
    element_t& operator[](int index) {return c.vec[index]->element;}
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

  // Provides an `Sequence` which can be used to efficiently iterate through the elements
  Sequence elements()
  {
    return *this;
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
