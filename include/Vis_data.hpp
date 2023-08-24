#ifndef HEXED_VIS_DATA_HPP_
#define HEXED_VIS_DATA_HPP_

#include "Element.hpp"
#include "Qpoint_func.hpp"
#include "Basis.hpp"

namespace hexed
{

/*!
 * Computes data to be visualized for a single element
 * (e.g. edge positions, non-conserved variables, values at uniformly-spaced sample points)
 * without knowing anything about the choice of visualization software.
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
  /*!
   * \param elem Element for which visualization is to be performed
   * \param func Computes variables to be visualized
   * \param basis Basis with which to perform interp/extrapolation
   * \param time Flow time used to compute output variables
   */
  Vis_data(Element& elem, const Qpoint_func& func, const Basis& basis, double time = 0.);
  /*! \brief interpolate function to `n_sample + 1` uniformly spaced points along element edges
   * \details layout: [number of edges in element][n_var (of Qpoint_func)][n_sample]
   */
  Eigen::VectorXd edges(int n_sample = 21);
  /*! \brief interpolate function to uniformly-spaced block of sample points `n_sample` on a side
   * \details layout: [n_var][n_sample]([n_sample]([n_sample]))
   */
  Eigen::VectorXd interior(int n_sample = 21);
  //! \brief interpolate function t a uniformly-spaced block of sample_points on a specified face
  Eigen::VectorXd face(int i_dim, bool is_positive, int n_sample = 21);
  //! \brief return function evaluated at quadratre points. \details layout: [n_var][n_qpoint]
  inline const Eigen::VectorXd& qpoints() {return vars;}

  //! \brief stores data representing a contour line/surface, to be converted to an otter curve/surface object
  struct Contour
  {
    Eigen::MatrixXd vert_ref_coords; //!< coordinates of contour vertices in reference coordinates
    Eigen::MatrixXd normals; //!< unit normal vectors to the contour surface in physical space (not reference) located at vertices
    Eigen::MatrixXi elem_vert_inds; //!< indices of contour elements (line segments/quads). layout: [i_element][i_vertex]
  };
  /*! \brief sample the function at a set of points given in reference coordinates
   * \param ref_coords: reference coordinates of sample points. layout: [n_dim][n_sample]
   * \returns values of visualization variables at sample points. layout: [n_var][n_sample]
   */
  Eigen::MatrixXd sample(Eigen::MatrixXd ref_coords);
  /*! \brief compute a contour line/surface where the `i_var`th variable is equal to `value`
   * \details the number of sample points in each direction is `2*n_div + 1`
   */
  Contour compute_contour(double value, int n_div = 10, int n_newton = 4, double tol = 1e-3);
};

}
#endif
