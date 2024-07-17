#ifndef HEXED_GEOM_EDGE_HPP_
#define HEXED_GEOM_EDGE_HPP_

#include <string>
#include "Array.hpp"

namespace hexed {

class Geom_edge {
  public:
  Geom_edge(Array<double>&& points);
  inline const Array<double> points() const {return _points();}
  inline const Array<double> arc_len() const {return _arc_len();}
  void visualize(std::string format, std::string name) const;

  private:
  Array<double> _points;
  int _n_points;
  Array<double> _arc_len;
};

}
#endif
