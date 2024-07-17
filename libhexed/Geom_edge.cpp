#include <hexed/Geom_edge.hpp>
#include <hexed/Visualizer.hpp>

namespace hexed {

Geom_edge::Geom_edge(Array<double>&& p)
: _points{std::move(p)}, _n_points{_points.shape()[0]}, _arc_len({_n_points})
{
  HEXED_ASSERT(_points.order() == 2, "point array must be order 2");
  HEXED_ASSERT(_points.shape()[1] == 3, "point array must have 3 columns");
  HEXED_ASSERT(_n_points >= 2, "point array must have at least 2 rows");
  _arc_len[0] = 0;
  for (int i_point = 1; i_point < _n_points; ++i_point) {
    _arc_len[i_point] = _arc_len[i_point - 1] + (_points(i_point) - _points(i_point - 1)).vector().norm();
  }
}

void Geom_edge::visualize(std::string format, std::string name) const {
  auto vis = Visualizer::create(format, 3, 1, name, std::vector<std::string>{}, 0., Visualizer::block);
  Array<double> transposed({3, _n_points});
  for (int i_point = 0; i_point < _n_points; ++i_point) {
    for (int i_dim = 0; i_dim < 3; ++i_dim) transposed(i_dim)[i_point] = _points(i_point)[i_dim];
  }
  vis->write_block(transposed, Array<double>({0, _n_points}));
}

Geom_edge::Node Geom_edge::nearest(Mat<3> to, double start, double stop) const {
  Node node{Mat<3>::Constant(std::nan("")), 0, std::nan("")};
  double dist = huge;
  for (int i_point = 0; i_point < _n_points; ++i_point) {
    double a = _arc_len[i_point];
    if (a >= start && a < stop) {
      Mat<3> p = _points(i_point).vector();
      double d = (p - to).norm();
      if (d < dist) {
        dist = d;
        node.pos = p;
        node.index = i_point;
        node.arc_len = a;
      }
    }
  }
  return node;
}

}
