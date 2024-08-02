#ifndef HEXED_CAD_GEOM_HPP_
#define HEXED_CAD_GEOM_HPP_

#include "math.hpp"

namespace hexed::cad_geom {

template <int n_param>
class Entity {
  public:
  virtual Mat<3> temp_point(Mat<n_param>) const = 0;
  Mat<3> point(Mat<n_param> params) const {return translate + transform*temp_point(params);}
  Mat<3, 3> transform = Mat<3, 3>::Identity();
  Mat<3> translate = Mat<3>::Zero();
};

class Circular_arc : public Entity<1> {
  public:
  Mat<3> temp_point(Mat<1>) const override;
  Mat<3> center;
  double radius;
  double start_angle;
  double end_angle;
};

}
#endif
