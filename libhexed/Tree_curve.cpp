#include <hexed/Tree_curve.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

Array<double>& check(Array<double>& nodes) {
  HEXED_ASSERT(nodes.shape()[0] == math::pow(2, math::log(2, nodes.shape()[0] - 1)) + 1,
               "number of `nodes` must be one greater than a power of 2");
  HEXED_ASSERT(nodes.shape()[1] == 3, "`nodes` must have 3 columns");
  return nodes;
}

Tree_curve::Tree_curve(Array<double>&& nodes, int skip)
: _nodes(std::move(check(nodes)))
, _arc_length({_nodes.shape()[0]})
, _skip{skip}
, _levels{std::max<Int>(1, math::log(2, _nodes.shape()[0] - 1) - _skip)}
, _segments({math::pow(2, _levels) - 1})
{
  _arc_length[0] = 0.;
  for (int i_node = 1; i_node < _arc_length.size(); ++i_node) {
    _arc_length[i_node] = _arc_length[i_node - 1] + (_nodes(i_node).vector() - _nodes(i_node - 1).vector()).norm();
  }
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
                           _segments(inds[0], inds[1]), n(), n_div*i_segment);
    }
  }
}

const Array<Tree_curve::Segment> Tree_curve::segments(int level) const {
  Int p = math::pow<Int>(2, level);
  return _segments(p - 1, 2*p - 1);
}

Tree_curve::Nearest_index Tree_curve::nearest_point(Mat<3> point, double max_dist, std::array<double, 2> bounds) const {
  Nearest_index nearest(-1, max_dist, -1.);
  _recursive_nearest(point, nearest, root(), bounds);
  return nearest;
}

void Tree_curve::_recursive_nearest(Mat<3> point, Nearest_index& nearest, const Segment& s,
                                    std::array<double, 2> bounds) const {
  if (s.segments.size()) {
    // call _recursive_nearest on both of the child segments if they are close enough
    // that they could possibly contain the nearest point
    double dist [2];
    for (int i_segment = 0; i_segment < 2; ++i_segment) {
      dist[i_segment] = (s.segments[i_segment].center - point).norm();
    }
    // do the closer segment first in hopes that we can find a point close enough to justify skipping the farther one
    for (bool i_segment : {dist[1] < dist[0], !(dist[1] < dist[0])}) {
      if (dist[i_segment] - s.segments[i_segment].radius < nearest.distance) {
        _recursive_nearest(point, nearest, s.segments[i_segment], bounds);
      }
    }
  } else {
    HEXED_ASSERT(s.nodes.shape()[0] >= 2, "`Segment` should have at least 2 nodes");
    for (Int i_node = 0; i_node < s.nodes.shape()[0] - 1; ++i_node) {
      Mat<3> diff = s.nodes(i_node + 1).vector() - s.nodes(i_node).vector();
      Mat<3> relative_point = point - s.nodes(i_node).vector();
      double interp = std::max(0., std::min(1., relative_point.dot(diff)/diff.squaredNorm()));
      double d = (relative_point - interp*diff).norm();
      double arc = (1 - interp)*_arc_length[s.nodes_start + i_node]
                       + interp*_arc_length[s.nodes_start + i_node + 1];
      if (d < nearest.distance && bounds[0] < arc && arc < bounds[1]) {
        nearest.interp_index = s.nodes_start + i_node + interp;
        nearest.index = std::round(nearest.interp_index);
        nearest.distance = d;
      }
    }
  }
}

//! \todo test this
Mat<3> Tree_curve::interp_point(double interp_index) const {
  Int index = std::max<Int>(0, std::min<Int>(_nodes.shape()[0] - 2, std::floor(interp_index)));
  double interp = interp_index - index;
  return (1 - interp)*_nodes(index).vector() + interp*_nodes(index + 1).vector();
}

Mat<3> Tree_curve::interp_point(Nearest_index near) const {
  HEXED_ASSERT(near.index >= 0, "`Nearest_index` indicates nonexistant node")
  return interp_point(near.interp_index);
}

std::vector<double> Tree_curve::intersections_2d(Mat<3> p0, Mat<3> p1) const {
  std::vector<double> sections;
  _recursive_intersections(p0, p1, sections, root());
  return sections;
}

void Tree_curve::_recursive_intersections(Mat<3> p0, Mat<3> p1, std::vector<double>& vec, const Segment& seg) const {
  auto seq = Eigen::seqN(0, 2);
  Mat<2> diff = p1(seq) - p0(seq);
  Mat<2> center_diff = seg.center(seq) - p0(seq);
  double distance_sq = (center_diff - center_diff.dot(diff)/diff.squaredNorm()*diff).squaredNorm();
  if (distance_sq <= seg.radius*seg.radius) {
    if (seg.segments.size()) {
      for (int i_segment = 0; i_segment < 2; ++i_segment) {
        _recursive_intersections(p0, p1, vec, seg.segments[i_segment]);
      }
    } else {
      for (int i_node = 0; i_node < seg.nodes.shape()[0] - 1; ++i_node) {
        Mat<2, 2> lhs;
        lhs(all, 0) = diff;
        lhs(all, 1) = (seg.nodes(i_node)(0, 2) - seg.nodes(i_node + 1)(0, 2)).vector();
        Mat<2> rhs = seg.nodes(i_node)(0, 2).vector() - p0(seq);
        Mat<2> soln = lhs.householderQr().solve(rhs);
        if (0 <= soln(1) && soln(1) <= 1) vec.push_back(soln(0));
      }
    }
  }
}

}
