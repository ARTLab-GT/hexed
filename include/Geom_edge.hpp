#ifndef HEXED_GEOM_EDGE_HPP_
#define HEXED_GEOM_EDGE_HPP_

#include "Array.hpp"

namespace hexed {

class Geom_edge {
  public:
  Geom_edge(Array<double>&& points);
  inline const Array<double> points() const {return _points();}
  inline const Array<double> arc_len() const {return _arc_len();}

  private:
  Array<double> _points;
  Array<double> _arc_len;
};

}
#endif
