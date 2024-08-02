#include <hexed/cad_geom.hpp>

namespace hexed::cad_geom {

Mat<3> Circular_arc::temp_point(Mat<1> params) const {
  double angle = start_angle + params(0)*(start_angle - end_angle);
  return center + radius*Mat<3>{std::cos(angle), std::sin(angle), 0.};
}

}
