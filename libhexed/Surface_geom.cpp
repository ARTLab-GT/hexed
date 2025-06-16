#include <hexed/Surface_geom.hpp>

namespace hexed {

Compound_edge::Compound_edge(std::vector<std::shared_ptr<Geom_edge>> edges, std::vector<bool> reverse)
: _edges(edges)
, _reverse(reverse)
{
  _reverse.resize(_edges.size(), false);
  double arc = 0;
  for (Int i_edge = 0; i_edge < (Int)_edges.size(); ++i_edge) {
    _arc_length_start.push_back(arc);
    arc += _edges[i_edge]->arc_length(1.);
  }
  _arc_length_start.push_back(arc);
}

Mat<3> Compound_edge::point(double param) const {
  auto ind = _global_index(param);
  return _edges[ind.first]->point(ind.second);
}

double Compound_edge::arc_length(double param) const {
  auto ind = _global_index(param);
  bool rev = _reverse[ind.first];
  return _arc_length_start[ind.first + rev] + math::sign(!rev)*_edges[ind.first]->arc_length(ind.second);
}

Mat<3> Compound_edge::tangent_average(double param) const {
  auto ind = _global_index(param);
  return _edges[ind.first]->tangent_average(ind.second);
}

double Compound_edge::tangent_radius(double param) const {
  auto ind = _global_index(param);
  return _edges[ind.first]->tangent_radius(ind.second);
}

double Compound_edge::arg_nearest_point(Mat<3> point) const {
  double dist_sq = huge;
  double param = 0;
  for (Int i_edge = 0; i_edge < (Int)_edges.size(); ++i_edge) {
    double p = _edges[i_edge]->arg_nearest_point(point);
    double dsq = (point - _edges[i_edge]->point(p)).squaredNorm();
    if (dsq < dist_sq) {
      dist_sq = dsq;
      if (_reverse[i_edge]) p = 1 - p;
      param = (p + i_edge)/_edges.size();
    }
  }
  return param;
}

std::pair<Int, double> Compound_edge::_global_index(double param) const {
  Int n = _edges.size();
  param *= n;
  Int edge = std::max<Int>(0, std::min<Int>(n - 1, floor(param)));
  param -= edge;
  if (_reverse[edge]) param = 1 - param;
  return {edge, param};
}

Compound_geom::Compound_geom(std::vector<Surface_geom*> geoms)
: components(geoms.begin(), geoms.end())
{}

Nearest_point<dyn> Compound_geom::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  Nearest_point nearest(point, max_distance);
  for (auto& comp : components) nearest.merge(comp->nearest_point(point, max_distance, distance_guess));
  return nearest;
}

std::vector<double> Compound_geom::intersections(Mat<> point0, Mat<> point1, bool high_precision) {
  std::vector<double> inters;
  for (auto& comp : components) {
    auto comp_inters = comp->intersections(point0, point1, high_precision);
    inters.insert(inters.end(), comp_inters.begin(), comp_inters.end());
  }
  return inters;
}

next::Sequence<const Geom_edge&> Compound_geom::edges() {
  next::Sequence<const Geom_edge&> e;
  for (auto& comp : components) e = e + comp->edges();
  return e;
}

next::Sequence<Mat<3>> Compound_geom::points() {
  next::Sequence<Mat<3>> p;
  for (auto& comp : components) p = p + comp->points();
  return p;
}

Hypersphere::Hypersphere(Mat<> center, double radius)
: c{center}, r{radius}
{}

Nearest_point<dyn> Hypersphere::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  point = resize(point, c.size());
  Nearest_point<dyn> nearest(point, max_distance);
  nearest.merge(c + r*(point - c).normalized());
  return nearest;
}

std::vector<double> Hypersphere::intersections(Mat<> point0, Mat<> point1, bool) {
  Mat<> start = point0 - c;
  Mat<> diff = point1 - point0;
  // a t^2 + b t + c = 0
  double a = diff.squaredNorm();
  double b = 2*diff.dot(start);
  double c = start.squaredNorm() - r*r;
  double discr = b*b - 4*a*c;
  if (discr < 0) return {};
  else if (discr == 0) return {-b/2/a};
  else {
    double root = std::sqrt(discr);
    return {(-b - root)/2/a, (-b + root)/2/a};
  }
}

}
