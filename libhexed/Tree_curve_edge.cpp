#include <hexed/Tree_curve_edge.hpp>

namespace hexed {

Tree_curve_edge::Tree_curve_edge(Array<double>&& nodes, Array<double>&& average, Array<double>&& radius,
                                 int skip_levels)
: _curve(std::move(nodes), skip_levels)
, _average{std::move(average)}
, _radius{std::move(radius)}
{}

Mat<3> Tree_curve_edge::point(double param) {
  return _curve.interp_point(param*(_curve.n_points() - 1));
}

double Tree_curve_edge::arc_length(double param) {
  return _curve.arc_length().flat_interp(param*(_curve.n_points() - 1));
}

Mat<3> Tree_curve_edge::tangent_average(double param) {
  return _average.interp(param*(_curve.n_points() - 1)).vector();
}

double Tree_curve_edge::tangent_radius(double param) {
  return _radius.flat_interp(param*(_curve.n_points() - 1));
}

double Tree_curve_edge::arg_nearest_point(Mat<3> p) {
  auto index = _curve.nearest_point(p);
  if (index.index >= 0) return index.interp_index;
  return 0.;
}

}
