#ifndef CARTDG_DEFORMED_ELEMENT_HPP_
#define CARTDG_DEFORMED_ELEMENT_HPP_

#include "Element.hpp"

namespace cartdg
{

/*
 * Represents an Element which is not a perfect square/cube. Note: Jacobian matrix
 * is nontrivial.
 */
class Deformed_element : public Element
{
  int n_qpoint;
  Eigen::VectorXd jac;
  Eigen::VectorXd node_adj;

  public:
  bool degenerate = 0;

  /*
   * `pos` sets the position of vertex 0 in increments of `mesh_size`. Unspecified
   * components set to 0. `mesh_size` sets distance between vertices. Vertex positions
   * may be modified later.
   */
  Deformed_element(Storage_params, std::vector<int> pos={}, double mesh_size=1.);
  virtual std::vector<double> position(const Basis&, int i_qpoint);
  // sets the Jacobian based on the current vertex locations and face node adjustments
  void set_jacobian(const Basis& basis);
  // Pointer to jacobian data. Use this for performance-cricial applications.
  double* jacobian(); // Layout: [i_dim][j_dim][i_qpoint]
  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  // No simple way to explain what this represents.
  double* node_adjustments(); // Layout: [i_dim][is_positive][i_face_qpoint]
};

}
#endif
