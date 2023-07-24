#include <Surface_geom.hpp>

namespace hexed
{

Compound_geom::Compound_geom(std::vector<Surface_geom*> geoms)
: components(geoms.begin(), geoms.end())
{}

math::Nearest_point<dyn> Compound_geom::nearest_point(Mat<> point, double max_distance, double distance_guess)
{
  math::Nearest_point nearest(point);
  for (auto& comp : components) nearest.merge(comp->nearest_point(point));
  return nearest;
}

std::vector<double> Compound_geom::intersections(Mat<> point0, Mat<> point1)
{
  std::vector<double> inters;
  for (auto& comp : components) {
    auto comp_inters = comp->intersections(point0, point1);
    inters.insert(inters.end(), comp_inters.begin(), comp_inters.end());
  }
  return inters;
}

Hypersphere::Hypersphere(Mat<> center, double radius)
: c{center}, r{radius}
{}

math::Nearest_point<dyn> Hypersphere::nearest_point(Mat<> point, double max_distance, double distance_guess)
{
  return math::Nearest_point<dyn>(point, c + r*(point - c).normalized());
}

std::vector<double> Hypersphere::intersections(Mat<> point0, Mat<> point1)
{
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
