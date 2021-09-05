#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <Eigen/Dense>

#include <Storage_params.hpp>

namespace cartdg
{

class Element
{
  unsigned n_stage;
  unsigned n_dof;
  double* data;
  void allocate();

  public:
  Element(Storage_params);
  Element(const Element& other);
  ~Element();
  Element& operator=(const Element& other);
  double* stage(unsigned i_stage);
};

}

#endif
