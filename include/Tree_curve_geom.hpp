#ifndef HEXED_TREE_CURVE_GEOM_HPP_
#define HEXED_TREE_CURVE_GEOM_HPP_

#include "Surface_geom.hpp"
#include "Tree_curve.hpp"

namespace hexed {

class Tree_curve_geom : public Surface_geom {
  public:
  Tree_curve_geom(Array<double>&& nodes, int skip_levels = 0);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  //! \warning Only 2D! Ignores the any entries past the second and assumes the \f$ x_2 \f$ coordinate is zero.
  std::vector<double> intersections(Mat<> point0, Mat<> point1, bool high_prec = true) override;
  next::Sequence<Mat<3>> points() override;

  private:
  Mat<3> _points_get(Int ind);
  Int _points_size();
  Tree_curve _curve;
};

}
#endif
