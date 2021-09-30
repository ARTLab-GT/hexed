#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Storage_params.hpp"

namespace cartdg
{

class Element
{
  int n_stage;
  int n_dof;
  int n_vert;
  Eigen::VectorXd data;
  Eigen::VectorXd visc_storage;
  Eigen::VectorXd derivative_storage;

  protected:
  int n_dim;

  public:
  Element(Storage_params);
  double* stage(int i_stage);
  double* face();
  // following two functions are for convenience, not performance
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  double jacobian_determinant(int i_qpoint);
  double* viscosity();
  bool viscous();
  double* derivative();
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::array<Element*, 2> elem_con;
typedef std::vector<std::vector<elem_con>> elem_con_vec;

}

#endif
