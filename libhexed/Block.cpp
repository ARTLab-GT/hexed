#include <hexed/Block.hpp>

namespace hexed::next
{

Mat<3> Block::point(std::vector<int> node_coords) const
{
  #ifdef DEBUG
  HEXED_ASSERT(int(node_coords.size()) == n_dim, "wrong number of node coordinates");
  for (int coord : node_coords) {
    HEXED_ASSERT(0 <= coord && coord < row_size, format_str(1000, "node index %i is out of bounds", coord));
  }
  #endif
  return _point(node_coords);
}

void Vertex::eat(Vertex& other) {}
void Vertex::average(std::vector<Vertex*>) {}

void Vertex::pair(Mutual_ptr<Edge, Vertex>& ptr)
{
  _edges.emplace_back(this);
  _edges.back().pair(ptr);
}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis> basis)
: Block(1, basis->row_size), _rs{basis->row_size}, _verts{this, this}, _interior({_rs, 3}), _basis{basis}
{
  vertex0.pair(_verts[0]);
  vertex1.pair(_verts[1]);
  reset();
}

Mat<3> Edge::_point(std::vector<int> coords) const
{
  int coord = coords[0];
  if (coord ==       0) return _verts[0]->point({});
  if (coord == _rs - 1) return _verts[1]->point({});
  return _interior(coord - 1).vector();
}

void Edge::reset()
{
  for (int i = 1; i < row_size - 1; ++i) {
    double n = _basis->node(i);
    _interior(i - 1) = Array<double>(Mat<3>((1 - n)*_verts[0]->point({}) + n*_verts[1]->point({})));
  }
}

void Edge::glue(Edge& other, double start, double stop)
{
}

}
