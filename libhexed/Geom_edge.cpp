#include <hexed/Geom_edge.hpp>

namespace hexed {

Geom_edge::Geom_edge(Array<double>&& p)
: _points{p}, _arc_len({_points.shape()[0]})
{
  HEXED_ASSERT(_points.order() == 2, "point array must be order 2");
  HEXED_ASSERT(_points.shape()[1] == 3, "point array must have 3 columns");
  HEXED_ASSERT(_points.shape()[0] >= 2, "point array must have at least 2 rows");
  _arc_len[0] = 0;
  for (int i_point = 1; i_point < _points.shape()[0]; ++i_point) {
    _arc_len[i_point] = _arc_len[i_point - 1] + (_points(i_point) - _points(i_point - 1)).vector().norm();
  }
}

}
