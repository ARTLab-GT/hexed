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

  protected:
  int n_dim;

  public:
  Element(Storage_params);
  double* stage(int i_stage);
  // following two functions are for convenience, not performance
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  double jacobian_determinant(int i_qpoint);
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::array<Element*, 2> elem_con;
typedef std::vector<std::vector<elem_con>> elem_con_vec;

}

#endif
