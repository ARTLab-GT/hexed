#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <Storage_params.hpp>

namespace cartdg
{

class Element
{
  unsigned n_stage;
  unsigned n_dof;
  Eigen::VectorXd data;

  public:
  Element(Storage_params);
  double* stage(unsigned i_stage);
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;

}

#endif
