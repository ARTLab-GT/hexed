#include <Simplex_geom.hpp>

namespace hexed
{

template<>
void Simplex_geom<2>::merge(Nearest_point<2>& nearest, Mat<2, 2> sim, Mat<2> point)
{
  nearest.merge(math::proj_to_segment({sim(all, 0), sim(all, 1)}, point));
}

template<>
void Simplex_geom<3>::merge(Nearest_point<3>& nearest, Mat<3, 3> sim, Mat<3> point)
{
  // try projecting the point to the plane of the triangle
  Mat<3, 2> lhs;
  lhs(all, 0) = sim(all, 1) - sim(all, 0);
  lhs(all, 1) = sim(all, 2) - sim(all, 0);
  Mat<2> lstsq = (lhs.transpose()*lhs).inverse()*(lhs.transpose()*(point - sim(all, 0)));
  if (lstsq(0) >= 0 && lstsq(1) >= 0 && lstsq.sum() <= 1) {
    // if the projected point is inside the triangle, evaluate it as the potential nearest point
    nearest.merge(sim(all, 0) + lhs*lstsq);
  } else {
    // if the projected point is outside the triangle, fall back to finding the nearest point on all the edges of the triangle
    for (int i_edge = 0; i_edge < 3; ++i_edge) {
      nearest.merge(math::proj_to_segment({sim(all, i_edge), sim(all, (i_edge + 1)%3)}, point));
    }
  }
}

std::vector<Mat<2, 2>> segments(const Mat<dyn, dyn>& points)
{
  std::vector<Mat<2, 2>> sims;
  for (unsigned i = 0; i < points.cols() - 1; ++i) {
    Mat<2, 2> sim;
    sim(all, 0) = points(all, i);
    sim(all, 1) = points(all, i + 1);
    sims.push_back(sim);
  }
  return sims;
}

}
