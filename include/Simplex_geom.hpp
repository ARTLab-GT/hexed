#ifndef HEXED_SIMPLEX_GEOM_HPP_
#define HEXED_SIMPLEX_GEOM_HPP_

#include "Surface_geom.hpp"

// At time of writing, this file is not `#include`d in any `.cpp` files,
// so the following command is necessary to make Doxygen parse it.
//! \file Simplex_geom.hpp

namespace hexed
{

template <int n_dim>
class Simplex_geom : public Surface_geom
{
  std::vector<Mat<n_dim, n_dim>> elems;

  public:
  Simplex_geom(const std::vector<Mat<n_dim, n_dim>>&  elements) : elems{elements} {}
  Simplex_geom(      std::vector<Mat<n_dim, n_dim>>&& elements) : elems{elements} {}

  Mat<> nearest_point(Mat<> point) override;

  std::vector<double> intersections(Mat<> point0, Mat<> point1) override
  {
    std::vector<double> inters;
    Mat<n_dim> diff = point1 - point0;
    for (Mat<n_dim, n_dim> elem : elems) {
      Mat<n_dim, n_dim> lhs;
      lhs(all, 0) = -diff;
      for (int col = 1; col < n_dim; ++col) lhs(all, col) = elem(all, col) - elem(all, 0);
      Mat<n_dim> soln = lhs.householderQr().solve(point0 - elem(all, 0));
      Eigen::Array<double, n_dim - 1, 1> arr = soln(Eigen::seqN(1, n_dim - 1)).array();
      if ((arr >= 0.).all() && arr.sum() <= 1.) inters.push_back(soln(0));
    }
    return inters;
  }
};

template<>
Mat<> Simplex_geom<2>::nearest_point(Mat<> point)
{
  math::Nearest_point<2> nearest(point);
  for (Mat<2, 2> elem : elems) {
    nearest.merge(math::proj_to_segment({elem(all, 0), elem(all, 1)}, point));
  }
  return nearest.point();
}

template <>
Mat<> Simplex_geom<3>::nearest_point(Mat<> point)
{
  math::Nearest_point<3> nearest(point);
  for (Mat<3, 3> elem : elems) {
    // try projecting the point to the plane of the triangle
    Mat<3, 2> lhs;
    lhs(all, 0) = elem(all, 1) - elem(all, 0);
    lhs(all, 1) = elem(all, 2) - elem(all, 0);
    Mat<2> lstsq = lhs.householderQr().solve(point - elem(all, 0));
    if (lstsq(0) >= 0 && lstsq(1) >= 0 && lstsq.sum() <= 1) {
      // if the projected point is inside the triangle, evaluate it as the potential nearest point
      nearest.merge(elem(all, 0) + lhs*lstsq);
    } else {
      // if the projected point is outside the triangle, fall back to finding the nearest point on all the edges of the triangle
      for (int i_edge = 0; i_edge < 3; ++i_edge) {
        nearest.merge(math::proj_to_segment({elem(all, i_edge), elem(all, (i_edge + 1)%3)}, point));
      }
    }
  }
  return nearest.point();
}

}
#endif
