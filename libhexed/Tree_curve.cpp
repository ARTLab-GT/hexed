#include <hexed/Tree_curve.hpp>

namespace hexed {

Array<double>& check(Array<double>& nodes) {
  HEXED_ASSERT(nodes.shape()[0] == math::pow(2, math::log(2, nodes.shape()[0] - 1)) + 1,
               "number of `nodes` must be one greater than a power of 2");
  HEXED_ASSERT(nodes.shape()[1] == 3, "`nodes` must have 3 columns");
  return nodes;
}

Tree_curve::Tree_curve(Array<double> nodes, int skip)
: _nodes(std::move(check(nodes)))
, _skip{skip}
, _levels{math::log(2, _nodes.shape()[0] - 1) - skip}
, _segments({math::pow(2, _levels) - 1})
{
  for (int level = _levels - 1; level >= 0; --level) {
    Int n_segs = math::pow(2, level);
    Int n_div = (_nodes.shape()[0] - 1)/n_segs;
    for (Int i_segment = 0; i_segment < n_segs; ++i_segment) {
      Array<double> n(_nodes(n_div*i_segment, n_div*(i_segment + 1) + 1));
      Mat<3> center = Mat<3>::Zero();
      for (Int i_node = 0; i_node < n_div + 1; ++i_node) center += n(i_node).vector();
      center /= n_div + 1;
      double rsq = 0;
      for (Int i_node = 0; i_node < n_div + 1; ++i_node) {
        rsq = std::max(rsq, (center - n(i_node).vector()).squaredNorm());
      }
      Int inds [2] {2*(n_segs + i_segment) - 1, 2*(n_segs + i_segment) + 1};
      if (level == _levels - 1) inds[0] = inds[1] = 0;
      _segments.initialize(n_segs - 1 + i_segment, center, std::sqrt(rsq),
                           _segments(inds[0], inds[1]), n, n_div*i_segment);
    }
  }
}

const Array<Tree_curve::Segment> Tree_curve::segments(int level) const {
  Int p = math::pow<Int>(2, level);
  return _segments(p - 1, 2*p - 1);
}

Nearest_point<3> Tree_curve::nearest_point(Mat<3> point, double max_dist) const {
  Nearest_point<3> nearest(point, max_dist);
  _recursive_nearest(nearest, root());
  return nearest;
}

void Tree_curve::_recursive_nearest(Nearest_point<3>& p, const Segment& s) const {
  if (math::pow(std::max((s.center - p.reference()).norm() - s.radius, 0.), 2) <= p.dist_squared()) {
    if (s.segments.size()) {
      for (int i_segment = 0; i_segment < 2; ++i_segment) _recursive_nearest(p, s.segments[i_segment]);
    } else {
      for (Int i_node = 0; i_node < s.nodes.shape()[0]; ++i_node) {
        p.merge(s.nodes(i_node).vector());
      }
    }
  }
}

}
