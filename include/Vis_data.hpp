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
  Array<double> interior(Int n_sample) const;
  /*! \brief interpolate function to `n_sample + 1` uniformly spaced points along element edges
   * \details layout: [number of edges in element][n_var (of Qpoint_func)][n_sample]
   */
  Array<double> edges(Int n_sample) const;

  //! \brief stores data representing a contour line/surface
  struct Contour {
    Array<double> vert_ref_coords; //!< coordinates of contour vertices in reference coordinates. layout: [n_dim][n_vertex]
    Array<Int> elem_vert_inds; //!< indices of contour elements (line segments/quads). layout: [i_element][math::pow(2, n_dim - 1)]
  };
  /*! \brief compute a contour line/surface where the `i_var`th variable is equal to `value`
   * \details the number of sample points in each direction is `2*n_div + 1`
   */
  Contour compute_contour(int i_var, double value, Int n_div, int n_newton = 4, double tol = 1e-3);

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
