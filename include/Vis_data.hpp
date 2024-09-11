#ifndef HEXED_VIS_DATA_HPP_
#define HEXED_VIS_DATA_HPP_

#include "Basis.hpp"
#include "Array.hpp"

namespace hexed {

/*! \brief Computes data to be visualized for a single element
 * \details (e.g. edge positions, non-conserved variables, values at uniformly-spaced sample points)
 * without knowing anything about the choice of visualization software.
 */
class Vis_data {
  public:
  /*!
   * \param qpoint_data Data to be visualized, at the quadrature points.
   * \param basis Basis with which to perform interp/extrapolation
   */
  Vis_data(Array<double> qpoint_data, const Basis& basis);
  /*! \brief Samples the quadrature point data at a set of arbitrary reference coordinates.
   * \param coords Array of reference coordinates to sample at. Layout: [n_dim][n_sample_point]
   * \return Array of data sampled at the specified points. Layout: [n_var][n_sample_point]
   */
  Array<double> sample(Array<double> coords) const;
  #if 0
  /*! \brief interpolate function to `n_sample + 1` uniformly spaced points along element edges
   * \details layout: [number of edges in element][n_var (of Qpoint_func)][n_sample]
   */
  Eigen::VectorXd edges(int n_sample = 21);
  /*! \brief interpolate function to a uniformly-spaced block of sample points `n_sample` on a side
   * \details layout: [n_var][n_sample]([n_sample]([n_sample]))
   */
  Eigen::VectorXd interior(int n_sample = 21);
  //! \brief interpolate function to a uniformly-spaced block of sample points on a specified face
  Eigen::VectorXd face(int i_dim, bool is_positive, int n_sample = 21);
  //! \brief return function evaluated at quadrature points. \details layout: [n_var][n_qpoint]
  inline const Eigen::VectorXd& qpoints() {return vars;}

  //! \brief stores data representing a contour line/surface
  struct Contour {
    Eigen::MatrixXd vert_ref_coords; //!< coordinates of contour vertices in reference coordinates. layout: [n_vertex][n_dim]
    //! unit normal vectors to the contour surface in physical space (not reference) located at vertices. layout: [n_vertex][n_dim]
    Eigen::MatrixXd normals;
    Eigen::MatrixXi elem_vert_inds; //!< indices of contour elements (line segments/quads). layout: [i_element][math::pow(2, n_dim - 1)]
  };
  /*! \brief compute a contour line/surface where the `i_var`th variable is equal to `value`
   * \details the number of sample points in each direction is `2*n_div + 1`
   */
  Contour compute_contour(double value, int n_div = 10, int n_newton = 4, double tol = 1e-3);
  #endif

  private:
  Array<double> _sample(Array<double> data, Array<double> coords) const;
  Array<double> _data;
  Int _n_dim;
  Int _n_edge;
  Int _row_size;
  Int _n_qpoint;
  Int _n_var;
  const Basis& _basis;
};

}
#endif
