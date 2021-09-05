#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <Eigen/Dense>

#include <Storage_params.hpp>

namespace cartdg
{

class Element
{
  double* data;
  unsigned n_dof;

  public:
  Element(Storage_params);
  ~Element();
  // FIXME: copy and move
  double* stage(unsigned i_stage);
};

}

#endif
