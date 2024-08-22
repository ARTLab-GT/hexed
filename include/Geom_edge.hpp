#ifndef HEXED_GEOM_EDGE_HPP_
#define HEXED_GEOM_EDGE_HPP_

#include <string>
#include "Array.hpp"
#include "Block.hpp"

namespace hexed {

class Geom_edge {
  public:
  struct Node {
    Mat<3> pos;
    Int index;
    double arc_len;
  };

  std::vector<Mortal_ptr<next::Edge>> matched_edges;
  std::vector<Mortal_ptr<next::Vertex>> matched_vertices;

  Geom_edge(Array<double> points);
  inline const Array<double> points() const {return _points();}
  inline const Array<double> arc_len() const {return _arc_len();}
  inline Int n_points() const {return _n_points;}
  inline double len() const {return _arc_len[_arc_len.size()];}
  void visualize(std::string format, std::string name) const;
  Node nearest(Mat<3> to, double start = 0, double stop = huge) const;

  private:
  Array<double> _points;
  Int _n_points;
  Array<double> _arc_len;
};

}
#endif
