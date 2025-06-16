#include <hexed/Tree_curve_geom.hpp>

namespace hexed {

Tree_curve_geom::Tree_curve_geom(Array<double>&& nodes, int skip_levels)
: _curve(std::move(nodes), skip_levels)
{}

Nearest_point<dyn> Tree_curve_geom::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  int nd = std::min<int>(point.size(), 3);
  Mat<3> p = Mat<3>::Zero();
  for (int i = 0; i < nd; ++i) p(i) = point(i);
  Nearest_point<dyn> nearest(point, max_distance);
  auto index = _curve.nearest_point(p, max_distance);
  if (index.index >= 0) nearest.merge(resize(_curve.interp_point(index), point.size()));
  return nearest;
}

std::vector<double> Tree_curve_geom::intersections(Mat<> point0, Mat<> point1, bool high_prec) {
  Mat<3> p0 = Mat<3>::Zero();
  Mat<3> p1 = Mat<3>::Zero();
  for (int i = 0; i < std::min<int>(3, std::min(point0.size(), point1.size())); ++i) {
    p0(i) = point0(i);
    p1(i) = point1(i);
  }
  return _curve.intersections_2d(p0, p1);
}

next::Sequence<Mat<3>> Tree_curve_geom::points() {
  return {
    [this](Int ind)->Mat<3>{return _curve.nodes()(ind*(_curve.n_points() - 1)).vector();},
    []()->Int{return 2;}
  };
}

}
