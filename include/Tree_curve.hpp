#ifndef HEXED_TREE_CURVE_HPP_
#define HEXED_TREE_CURVE_HPP_

#include "Array.hpp"
#include "Nearest_point.hpp"

namespace hexed {

class Tree_curve {
  public:
  struct Segment {
    Mat<3> center;
    double radius;
    const Array<Segment> segments;
    const Array<double> nodes;
    Int nodes_start;
  };
  struct Nearest_index {
    Int index;
    double distance;
  };
  Tree_curve(Array<double>&& nodes, int skip_levels = 0);
  inline int skip_levels() const {return _skip;}
  inline const Array<double> nodes() const {return _nodes();}
  inline const Array<double> arc_length() const {return _arc_length();}
  const Array<Segment> segments(int level) const;
  inline const Segment& root() const {return _segments[0];}
  Nearest_index nearest_point(Mat<3> point, double max_dist = std::sqrt(huge),
                              std::array<double, 2> arc_len_bounds = {-huge, huge}) const;

  private:
  void _recursive_nearest(Mat<3> point, Nearest_index&, const Segment&, std::array<double, 2> arc_len_bounds) const;
  Array<double> _nodes;
  Array<double> _arc_length;
  int _skip;
  Int _levels;
  Array<Segment> _segments;
};

}
#endif
