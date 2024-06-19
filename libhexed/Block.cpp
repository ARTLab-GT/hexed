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

void Vertex::eat(Vertex& other)
{
  if (&other == this) return;
  other._alive = false;
  pos = (_mass*pos + other._mass*other.pos)/(_mass + other._mass);
  _mass += other._mass;
  other.purge();
  for (auto& edge : other._edges) pair(*edge.partner());
  other._mass = 0;
}

void Vertex::pair(Mutual_ptr<Edge, Vertex>& ptr)
{
  _edges.emplace_back(this);
  _edges.back().pair(ptr);
}

void Vertex::purge()
{
  std::erase(_edges, false);
}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis> basis) :
  Block(1, basis->row_size),
  _rs{basis->row_size},
  _verts{this, this},
  _interior({_rs, 3}),
  _basis{basis},
  _glued_to{this}
{
  vertex0.pair(_verts[0]);
  vertex1.pair(_verts[1]);
  reset();
}

Mat<3> Edge::_point(std::vector<int> coords) const
{
  int coord = coords[0];
  if (glued()) {
    if (_half == no) return _glued_to->_point(coords);
    else {
      Mat<3, dyn> pts(3, _rs);
      for (int c = 0; c < _rs; ++c) pts(all, c) = _glued_to->point({c});
      return pts*_basis->restrict(_half)(coord, all).transpose();
    }
  }
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

const int Edge::no = -1;

void Edge::glue(Edge& other, int half)
{
  other._glued.emplace_back(&other);
  _glued_to.pair(other._glued.back());
  _half = half;
}

void Edge::unglue()
{
  Edge& other = _glued_to.value();
  _glued_to.unpair();
  std::erase(other._glued, false);
}

bool Edge::glued() const
{
  return _glued_to;
}

Mat<3> Surface_face::_point(std::vector<int> coords) const
{
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    if (coords[i_dim] ==       0) return _edges[2*i_dim    ].point({coords[!i_dim]});
    if (coords[i_dim] == _rs - 1) return _edges[2*i_dim + 1].point({coords[!i_dim]});
  }
  return Mat<3>::Zero();
}

Surface_face::Surface_face(std::array<Vertex*, 4> verts, std::shared_ptr<Basis> basis)
: Block(2, basis->row_size), _rs{basis->row_size}, _interior({_rs*_rs, 3})
{
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _edges.emplace_back(*verts[(2 - i_dim)*sign], *verts[(2 - i_dim)*sign + 1 + i_dim], basis);
    }
  }
}

void Surface_face::reset()
{
}

}
