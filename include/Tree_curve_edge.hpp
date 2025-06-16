#ifndef HEXED_TREE_CURVE_EDGE_HPP_
#define HEXED_TREE_CURVE_EDGE_HPP_

#include "Tree_curve.hpp"
#include "Surface_geom.hpp"

namespace hexed {

class Tree_curve_edge : public Geom_edge {
  public:
  Tree_curve_edge(Array<double>&& nodes, Array<double>&& tangent_average, Array<double>&& tangent_radius,
                  int skip_levels = 0);
  Mat<3> point(double) const override;
  double arc_length(double) const override;
  Mat<3> tangent_average(double) const override;
  double tangent_radius(double) const override;
  double arg_nearest_point(Mat<3>) const override;
  private:
  Tree_curve _curve;
  Array<double> _average;
  Array<double> _radius;
};

}
#endif
