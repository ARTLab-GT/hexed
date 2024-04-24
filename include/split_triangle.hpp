#ifndef HEXED_SPLIT_TRIANGLE_HPP_
#define HEXED_SPLIT_TRIANGLE_HPP_

#include "math.hpp"

namespace hexed
{

//! \brief computes the area of a triangle
template <int n_dim>
double triangle_area(Mat<n_dim, 3> triangle)
{
  Mat<n_dim> base = triangle(all, 1) - triangle(all, 0);
  Mat<n_dim> side = triangle(all, 2) - triangle(all, 0);
  double bsq = base.squaredNorm();
  return 0.5*sqrt(bsq*(side - side.dot(base)/bsq*base).squaredNorm());
}

/*! \brief splits a triangle into 4 congruent parts
 * \param verts each column is one vertex of the triangle
 * \param params if provided, each column is the parametric coordinates of one vertex with respect to some analytic surface
 */
std::vector<std::pair<Mat<3, 3>, Mat<2, 3>>> split_triangle(Mat<3, 3> verts, Mat<2, 3> params = Mat<2, 3>::Zero());

}
#endif
