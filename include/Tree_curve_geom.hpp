#ifndef HEXED_TREE_CURVE_GEOM_HPP_
#define HEXED_TREE_CURVE_GEOM_HPP_

#include "Surface_geom.hpp"
#include "Tree_curve.hpp"

namespace hexed {

class Tree_curve_geom : public Surface_geom {
  public:
  Tree_curve_geom(Array<double>&& nodes, int skip_levels = 0);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
  next::Sequence<Mat<3>> points() override;

  private:
  Tree_curve _curve;
};

}
#endif
