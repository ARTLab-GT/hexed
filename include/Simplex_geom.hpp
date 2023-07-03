#ifndef HEXED_SIMPLEX_GEOM_HPP_
#define HEXED_SIMPLEX_GEOM_HPP_

#include "Surface_geom.hpp"

namespace hexed
{

template <int n_dim>
class Simplex_geom : public Surface_geom
{
  class Nearest_candidate
  {
    public:
    Mat<n_dim> point;
    double dist_sq;
    bool operator<(const Nearest_candidate& other) const {return dist_sq < other.dist_sq;}
  };

  std::vector<Mat<n_dim, n_dim>> elems;
  void project_to_elem(Mat<n_dim, n_dim>, Mat<n_dim>, Nearest_candidate&);

  public:
  Simplex_geom(const std::vector<Mat<n_dim, n_dim>>&  elements) : elems{elements} {}
  Simplex_geom(      std::vector<Mat<n_dim, n_dim>>&& elements) : elems{elements} {}

  Mat<> nearest_point(Mat<> point) override
  {
    Mat<n_dim> p = point;
    Nearest_candidate nearest{p, std::numeric_limits<double>::max()};
    for (Mat<n_dim, n_dim> elem : elems) project_to_elem(elem, p, nearest);
    return nearest.point;
  }

  std::vector<double> intersections(Mat<> point0, Mat<> point1) override
  {
    return {};
  }
};

template <>
void Simplex_geom<2>::project_to_elem(Mat<2, 2> elem, Mat<2> point, Nearest_candidate& nearest)
{
  Mat<2> proj = math::proj_to_segment({elem(all, 0), elem(all, 1)}, point);
  nearest = std::min(nearest, Nearest_candidate{proj, (proj - point).squaredNorm()});
}

template <>
void Simplex_geom<3>::project_to_elem(Mat<3, 3> elem, Mat<3> point, Nearest_candidate& nearest)
{
  Mat<3, 2> lhs;
  lhs(all, 0) = elem(all, 1) - elem(all, 0);
  lhs(all, 1) = elem(all, 2) - elem(all, 0);
  Mat<2> lstsq = lhs.householderQr().solve(point - elem(all, 0));
  if (lstsq(0) >= 0 && lstsq(1) >= 0 && lstsq.sum() <= 1) {
    Mat<3> proj = elem(all, 0) + lhs*lstsq;
    nearest = std::min(nearest, Nearest_candidate{proj, (proj - point).squaredNorm()});
  } else {
    for (int i_edge = 0; i_edge < 3; ++i_edge) {
      Mat<3> proj = math::proj_to_segment({elem(all, i_edge), elem(all, (i_edge + 1)%3)}, point);
      nearest = std::min(nearest, Nearest_candidate{proj, (proj - point).squaredNorm()});
    }
  }
}

}
#endif
