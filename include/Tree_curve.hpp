#ifndef HEXED_TREE_CURVE_HPP_
#define HEXED_TREE_CURVE_HPP_

#include "Array.hpp"

namespace hexed {

class Tree_curve {
  public:
  struct Segment {
    Mat<3> center;
    double radius;
    const Array<Segment> segments;
    const Array<double> nodes;
  };
  Tree_curve(Array<double> nodes, int skip_levels = 0);
  inline int skip_levels() const {return _skip;}
  inline const Array<double> nodes() const {return _nodes();}
  const Array<Segment> segments(int level) const;
  inline const Segment& root() const {return _segments[0];}

  private:
  Array<double> _nodes;
  int _skip;
  Int _levels;
  Array<Segment> _segments;
};

}
#endif
