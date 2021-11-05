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
 * This class represents a Cartesian (i.e., regular) element. See also derived class
 * `Deformed_element`.
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
  double* stage(int i_stage); // Layout: [i_var][i_qpoint]
  // Pointer state data at faces. Must be populated by user
  double* face(); // Layout: [i_dim][is_positive][i_var][i_face_qpoint]

  /*
   * Following two functions compute the Jacobian. I.e., derivative of `i_dim`th
   * physical coordinate wrt `j_dim`th reference coordinate. Trivial for this
   * class, may be non-trivial for derived (see `Deformed_element`). For
   * convenience, not performance.
   */
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint); // identity matrix
  double jacobian_determinant(int i_qpoint); // returns 1.

  double* viscosity(); // Artificial viscosity coefficient at corners.
  bool viscous(); // Should artificial viscosity be applied in this element?
  // Pointer to storage for derivative. Must be populated by user.
  double* derivative(); // Layout: [i_qpoint]
};

typedef std::vector<std::unique_ptr<Element>> elem_vec;
typedef std::array<Element*, 2> elem_con;
typedef std::vector<std::vector<elem_con>> elem_con_vec;

}

#endif
