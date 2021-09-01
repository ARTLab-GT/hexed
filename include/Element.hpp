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
  using ref_t = Eigen::VectorBlock<vec_t>;
  ref_t stage_block(int i_stage);
};

}

#endif
