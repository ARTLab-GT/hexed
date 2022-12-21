#ifndef HEXED_DEFORMED_ELEMENT_HPP_
#define HEXED_DEFORMED_ELEMENT_HPP_

#include "Element.hpp"

namespace hexed
{

/*
 * Represents an Element which is not a perfect axis-aligned square/cube.
 * Note: Jacobian matrix is nontrivial.
 */
class Deformed_element : public Element
{
  int n_qpoint;
  // jacobian data (first `n_dim*n_dim*n_qpoints` are `reference_level_normals`
  // and rest is `jacobian_determinant`)
  Eigen::VectorXd jac_dat;
  Eigen::VectorXd node_adj;
  double* f_nrml [6];

  public:
  bool degenerate = 0;
  static constexpr bool is_deformed = true;

  Deformed_element(Storage_params, std::vector<int> pos = {}, double mesh_size = 1., int ref_level = 0);
  virtual std::vector<double> position(const Basis&, int i_qpoint);
  // sets the Jacobian based on the current vertex locations and face node adjustments
  virtual void set_jacobian(const Basis& basis);

  /* the `j_dim`th component (in physical space) of the normal vector of the level surface
   * of the `i_dim`th reference coordinate which passes through the `i_qpoint`th quadrature point.
   * the magnitude of the normal vector is weighted by the surface area in physical space.
   * equivalent definition (note `i_dim`, `j_dim` transposed):
   * $$ J^{-1}_{j_dim, i_dim} |J| $$
   * where is jacobian of transformation from reference to physical coordinates at `i_qpoint`th quadrature point.
   * layout: [i_dim][j_dim][i_qpoint]
   */
  double* reference_level_normals();
  // jacobian determinant data
  double* jacobian_determinant();
  double*& face_normal(int i_face);

  virtual double jacobian(int i_dim, int j_dim, int i_qpoint);
  virtual double jacobian_determinant(int i_qpoint);
  // represents adjustments to face quadrature points to fit surfaces (details are complicated)
  virtual double* node_adjustments(); // Layout: [i_dim][is_positive][i_face_qpoint]
};

}
#endif
