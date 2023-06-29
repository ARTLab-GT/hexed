#include <Occt_geom.hpp>
#if HEXED_USE_OCCT

namespace hexed
{

Occt_geom::Occt_geom(TopoDS_Shape&&)
{}

Mat<> Occt_geom::nearest_point(Mat<> point)
{
  return Mat<>::Zero(3);
}

std::vector<double> Occt_geom::intersections(Mat<> point0, Mat<> point1)
{
  return {};
}

}
#endif
