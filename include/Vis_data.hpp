#ifndef HEXED_VIS_DATA_HPP_
#define HEXED_VIS_DATA_HPP_

#include "Element.hpp"
#include "Qpoint_func.hpp"
#include "Basis.hpp"

namespace hexed
{

/*
 * Computes data to be visualized (e.g. edge positions, non-conserved variables) without
 * knowing anything about the choice of visualization software.
 */
class Vis_data
{
  int n_dim;
  int n_edge;
  int row_size;
  int n_qpoint;
  int n_var;
  Element& el;
  const Basis& bas;
  Eigen::VectorXd vars;
  Eigen::MatrixXd sample_qpoint_data(Eigen::VectorXd qpoint_data, Eigen::MatrixXd ref_coords);

  public:
  Vis_data(Element&, const Qpoint_func&, const Basis&, double time = 0.);
  /*
   * interpolate function to `n_sample + 1` uniformly spaced points along element edges
   * layout: [number of edges in element][n_var (of Qpoint_func)][n_sample]
   */
  Eigen::VectorXd edges(int n_sample = 21);
  /*
   * interpolate function to uniformly-spaced block of sample points `n_sample` on a side
   * layout: [n_var][n_sample]([n_sample]([n_sample]))
   */
  Eigen::VectorXd interior(int n_sample = 21);
  // interpolate function t a uniformly-spaced block of sample_points on a specified face
  Eigen::VectorXd face(int i_dim, bool is_positive, int n_sample = 21);
  // return function evaluated at quadratre points. layout: [n_var][n_qpoint]
  inline const Eigen::VectorXd& qpoints() {return vars;}

  struct Contour
  {
    Eigen::MatrixXd vert_ref_coords;
    Eigen::MatrixXd normals; // normals in physical space (not reference)
    Eigen::MatrixXi elem_vert_inds;
  };
  // sample the function at a set of points given in reference coordinates
  // input layout: [n_dim][n_sample]
  // output layout: [n_var][n_sample]
  Eigen::MatrixXd sample(Eigen::MatrixXd ref_coords);
  // compute a contour line/surface where the `i_var`th variable is equal to `value`
  // the number of sample points in each direction is `2*n_div + 1`
  Contour compute_contour(double value, int n_div = 10, int n_newton = 4, double tol = 1e-3);
};

}
#endif
