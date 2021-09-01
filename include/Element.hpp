#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <Eigen/Dense>

#include <Storage_params.hpp>

namespace cartdg
{

class Element
{
  using vec_t = Eigen::VectorXd;
  vec_t data;
  int n_dof;

  public:
  Element(Storage_params);
  double* read_ptr();
  double* write_ptr();
};

}

#endif
