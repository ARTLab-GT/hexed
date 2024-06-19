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
  _interior({_rs - 2, 3}),
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
  return _interior(coords[0] - 1)(coords[0] - 1).vector();
}

Surface_face::Surface_face(std::array<Vertex*, 4> verts, std::shared_ptr<Basis> basis)
: Block(2, basis->row_size), _rs{basis->row_size}, _interior({_rs - 2, _rs - 2, 3}), _basis{basis}
{
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _edges.emplace_back(*verts[(2 - i_dim)*sign], *verts[(2 - i_dim)*sign + 1 + i_dim], basis);
    }
  }
  reset();
}

void Surface_face::reset()
{
  int sz = (_rs - 2)*(_rs - 2);
  Mat<dyn, dyn> dm = _basis->diff_mat();
  Mat<dyn, dyn> dmsq = dm*dm;
  std::cout << dmsq << "\n\n";
  std::cout << dmsq(Eigen::seqN(1, _rs - 2), Eigen::seqN(1, _rs - 2)).fullPivLu().rank() << "\n\n";
  Mat<dyn, dyn> lhs_mat(sz, sz);
  Mat<dyn, dyn> lhs_mat1(sz, sz);
  Mat<dyn, dyn> rhs_mat = Mat<dyn, dyn>::Zero(3, sz);
  Array<double> lhs(std::vector<int>(4, _rs - 2), lhs_mat.data());
  Array<double> lhs1(std::vector<int>(4, _rs - 2), lhs_mat1.data());
  Array<double> rhs({_rs - 2, _rs - 2, 3}, rhs_mat.data());
  for (int i_row = 0; i_row < _rs - 2; ++i_row) {
    for (int j_row = 0; j_row < _rs - 2; ++j_row) {
      for (int i_col = 0; i_col < _rs - 2; ++i_col) {
        for (int j_col = 0; j_col < _rs - 2; ++j_col) {
          lhs (i_row)(j_row)(i_col)[j_col] = dmsq(i_row + 1, i_col + 1);
          lhs1(i_row)(j_row)(i_col)[j_col] = dmsq(j_row + 1, j_col + 1);
        }
      }
      for (int sign = 0; sign < 2; ++sign) {
        int col = sign*(_rs - 1);
        int row [2] {i_row, j_row};
        for (int j_dim = 0; j_dim < 2; ++j_dim) {
          rhs(i_row)(j_row).vector() -= dmsq(row[j_dim], col)*edge(2*j_dim + sign).point({row[!j_dim]});
        }
      }
    }
  }
  Mat<dyn, dyn> lhs_mat2 = lhs_mat + lhs_mat1;
  std::cout << lhs_mat << "\n\n";
  std::cout << lhs_mat1 << "\n\n";
  std::cout << lhs_mat2 << "\n\n";
  std::cout << rhs_mat << "\n\n";
  Mat<dyn, dyn> rhs_t = rhs_mat.transpose();
  //Eigen::FullPivHouseholderQR<Mat<dyn, dyn>> fact(sz, sz);
  //fact.setThreshold(1e-16);
  //fact.compute(lhs_mat);
  {
  auto fact = lhs_mat.fullPivHouseholderQr();
  std::cout << fact.inverse() << "\n\n";
  }{
  auto fact = lhs_mat1.fullPivHouseholderQr();
  std::cout << fact.inverse() << "\n\n";
  }{
  auto fact = lhs_mat2.fullPivHouseholderQr();
  std::cout << fact.inverse() << "\n\n";
  }
  //HEXED_ASSERT(fact.rank() == sz, format_str(100, "LU factorization failed: rank = %i", fact.rank()));
  Mat<dyn, dyn> soln = lhs_mat.householderQr().solve(rhs_t);
  std::cout << soln << "\n\n";
  soln.transposeInPlace();
  soln.resize(3*sz, 1);
  _interior.vector() = soln;
}

}
