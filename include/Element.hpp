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
  int n_stage;
  int n_dof;
  Eigen::VectorXd data;

  public:
  Element(Storage_params);
  double* stage(int i_stage);
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::vector<std::vector<Element*>> elem_cons;

}

#endif
