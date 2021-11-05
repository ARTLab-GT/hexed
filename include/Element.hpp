#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "Storage_params.hpp"

namespace cartdg
{

/*
 * Stores data associated with one element. Container only --
 * does not have implementations of or information about the basis and algorithms.
 */
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
  // Pointer to state data for `i_stage`th Runge-Kutta stage.
  double* stage(int i_stage);
  double* face(); // Pointer state data for all faces. Must be populated by user
  // Following two functions compute the Jacobian of the transformation
  // from reference to physical coordinates. Trivial for this class, may be non-trivial
  // for derived (see `Deformed_element`). For convenience, not performance.
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint); // identity matrix
  double jacobian_determinant(int i_qpoint); // returns 1.
  double* viscosity(); // Artificial viscosity coefficient at corners.
  bool viscous(); // Should artificial viscosity be applied in this element?
  // Pointer to storage for derivative. Size: n_qpoint. Must be populated by user.
  double* derivative();
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::array<Element*, 2> elem_con;
typedef std::vector<std::vector<elem_con>> elem_con_vec;

}

#endif
