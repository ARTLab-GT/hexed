#include <hexed/Block.hpp>
#include <hexed/Visualizer.hpp>

namespace hexed::next
{

Array<double> Block::points() const
{
  std::vector<int> shape {3};
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) shape.push_back(row_size);
  Array<double> pts(shape);
  for (int i_point = 0; i_point < pts.stride(0); ++i_point) {
    std::vector<int> inds(n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) inds[i_dim] = (i_point/pts.stride(1 + i_dim))%row_size;
    auto pt = point(inds);
    for (int i_dim = 0; i_dim < 3; ++i_dim) pts(i_dim)[i_point] = pt(i_dim);
  }
  return pts;
}

void Block::visualize(std::string format, std::string file_name, const std::vector<Block*>& blocks, double time)
{
  int block_dim = blocks.empty() ? 1 : blocks[0]->n_dim;
  auto visualizer = Visualizer::create(format, 3, block_dim, file_name, {}, time, Visualizer::block);
  for (Block* block : blocks) {
    HEXED_ASSERT(block, "Null pointer passed to `Block::visualize`.");
    visualizer->write_block(block->points(), Array<double>({}));
  }
}

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

void Vertex::pair(Mutual_ptr<Mesh_element, Vertex>& ptr)
{
  _elems.emplace_back(this);
  _elems.back().pair(ptr);
}

void Vertex::purge()
{
  std::erase(_edges, false);
}

std::vector<int> interior_dims(int n_dim, int row_size)
{
  std::vector<int> dims(n_dim, row_size - 2);
  dims.push_back(3);
  return dims;
}

Boundary_interior::Boundary_interior(int n_dim, const Basis& b)
: Block(n_dim, b.row_size), _interior(interior_dims(n_dim, row_size)), basis{b}
{}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, const Basis& b) :
  Boundary_interior(1, b),
  _verts{this, this},
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
      Mat<3, dyn> pts(3, row_size);
      for (int c = 0; c < row_size; ++c) pts(all, c) = _glued_to->point({c});
      return pts*basis.restrict(_half)(coord, all).transpose();
    }
  }
  if (coord ==       0) return _verts[0]->point({});
  if (coord == row_size - 1) return _verts[1]->point({});
  return _interior(coord - 1).vector();
}

void Edge::reset()
{
  for (int i = 1; i < row_size - 1; ++i) {
    double n = basis.node(i);
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
    if (coords[i_dim] == row_size - 1) return _edges[2*i_dim + 1].point({coords[!i_dim]});
  }
  return _interior(coords[0] - 1)(coords[1] - 1).vector();
}

Surface_face::Surface_face(std::array<Vertex*, 4> verts, const Basis& b)
: Boundary_interior(2, b)
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
  int int_sz = (row_size - 2)*(row_size - 2);
  int tot_sz = row_size*row_size;
  Mat<dyn, dyn> dmsq = basis.diff_mat()*basis.diff_mat();
  Mat_rm<> lhs_mat = Mat_rm<>::Zero(tot_sz, int_sz);
  Mat_rm<> rhs_mat = Mat_rm<>::Zero(tot_sz, 3);
  Array<double> lhs({row_size, row_size, row_size - 2, row_size - 2}, lhs_mat.data());
  Array<double> rhs({row_size, row_size, 3}, rhs_mat.data());
  for (int i_row = 0; i_row < row_size; ++i_row) {
    for (int j_row = 1; j_row < row_size - 1; ++j_row) {
      for (int col = 1; col < row_size - 1; ++col) {
        lhs(i_row)(j_row)(col - 1)[j_row - 1] += dmsq(i_row, col);
        lhs(j_row)(i_row)(j_row - 1)[col - 1] += dmsq(i_row, col);
      }
      for (int col : {0, row_size - 1}) {
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(    bool(col)).point({j_row});
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(2 + bool(col)).point({j_row});
      }
    }
    for (int j_row : {0, row_size - 1}) {
      for (int col = 0; col < row_size; ++col) {
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(    bool(j_row)).point({col});
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(2 + bool(j_row)).point({col});
      }
    }
  }
  Mat_rm<> soln = lhs_mat.fullPivHouseholderQr().solve(rhs_mat);
  _interior = soln.data();
}

Mat<3> Mesh_element::_point(std::vector<int>) const
{
  return Mat<3>::Zero();
}

Mesh_element::Mesh_element(int nd, int rs)
: Block(nd, rs)
{
  for (int i_vert = 0; i_vert < math::pow(2, nd); ++i_vert) _verts.emplace_back(this);
}

const int Mesh_blocks::no_face = -1;

Mesh_blocks::Mesh_blocks(int nd, const Basis& b)
: n_dim{nd}, basis{b}
{}

Sequence<Vertex&> Mesh_blocks::interior_verts()
{
  return Sequence<Vertex&>::ptr_vector_view<std::unique_ptr<Vertex>&>(_interior_verts);
}

Sequence<Vertex&> Mesh_blocks::boundary_verts()
{
  return Sequence<Vertex&>::ptr_vector_view<std::unique_ptr<Vertex>&>(_boundary_verts);
}

std::unique_ptr<Mesh_element> Mesh_blocks::create_element(Mat<3> pos, double size, int boundary_face)
{
  std::unique_ptr<Mesh_element> ptr(new Mesh_element(n_dim, basis.row_size));
  int nv = math::pow(2, n_dim);
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    auto vec = &_interior_verts;
    if (boundary_face != no_face) {
      if ((i_vert/math::pow(2, n_dim - 1 - boundary_face/2))%2 == boundary_face%2) vec = &_boundary_verts;
    }
    vec->emplace_back(new Vertex(pos, basis.row_size));
    vec->back()->pair(ptr->_verts[i_vert]);
  }
  return ptr;
}

}
