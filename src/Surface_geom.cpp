#include <Surface_geom.hpp>

namespace hexed
{

Hypersphere::Hypersphere(Mat<> center, double radius)
: c{center}, r{radius}
{}

Mat<> Hypersphere::nearest_point(Mat<> point)
{
  return c + r*(point - c).normalized();
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
