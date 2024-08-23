#include <hexed/Tree_curve.hpp>

namespace hexed {

Tree_curve::Tree_curve(Array<double> nodes, int skip) : _nodes({}), _skip{skip}, _segments({}) {}

const Array<Tree_curve::Segment> Tree_curve::segments(int level) const {
  Int p = math::pow<Int>(2, level);
  //return _segments(p - 1, 2*p - 1);
  return _segments();
}

}
