#include <Surface_geom.hpp>

namespace hexed
{

Hypersphere::Hypersphere(Mat<> center, double radius)
: c{center}, r{radius}
{}

Mat<> Hypersphere::nearest_point(Mat<> point)
{
  return point;
}

std::vector<double> Hypersphere::intersections(Mat<> point0, Mat<> point1)
{
  return {};
}

}
