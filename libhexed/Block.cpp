#include <hexed/Block.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/vertex_inds.hpp>

namespace hexed::next {

Array<double> Block::points() const {
  std::vector<int> shape {3};
  for (int i_dim = 0; i_dim < _n_dim; ++i_dim) shape.push_back(_row_size);
  Array<double> pts(shape);
  for (int i_point = 0; i_point < pts.stride(0); ++i_point) {
    std::vector<int> inds(_n_dim);
    for (int i_dim = 0; i_dim < _n_dim; ++i_dim) inds[i_dim] = (i_point/pts.stride(1 + i_dim))%_row_size;
    auto pt = point(inds);
    for (int i_dim = 0; i_dim < 3; ++i_dim) pts(i_dim)[i_point] = pt(i_dim);
  }
  return pts;
}

void Block::visualize(std::string format, std::string file_name, const std::vector<Block*>& blocks, double time) {
  int block_dim = blocks.empty() ? 1 : blocks[0]->_n_dim;
  auto visualizer = Visualizer::create(format, 3, block_dim, file_name, {}, time, Visualizer::block);
  for (Block* block : blocks) {
    HEXED_ASSERT(block, "Null pointer passed to `Block::visualize`.");
    visualizer->write_block(block->points(), Array<double>({}));
  }
}

Mat<3> Block::point(std::vector<int> node_coords) const {
  #ifdef DEBUG
  HEXED_ASSERT(int(node_coords.size()) == _n_dim, "wrong number of node coordinates");
  for (int coord : node_coords) {
    HEXED_ASSERT(0 <= coord && coord < _row_size, format_str(1000, "node index %i is out of bounds", coord));
  }
  #endif
  return _point(node_coords);
}

Mat<3> Vertex::_point(std::vector<int>) const {
  if (!_glued_to) return pos;
  Array<double> points = _glued_to->points();
  Mat<3> p;
  Eigen::Map<const Mat<>> sample(_glued_coords.data(), _glued_coords.size());
  Mat<dyn, dyn> proj = _glued_to->basis().interpolate(sample);
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    Mat<> vec = points(i_dim).vector();
    for (int j_dim = _glued_coords.size() - 1; j_dim >= 0; --j_dim) {
      Mat<> new_vec = math::dimension_matvec(proj(j_dim, all), vec, j_dim);
      vec = new_vec;
    }
    p(i_dim) = vec(0);
  }
  return p;
}

Vertex::Vertex(Mat<3> pos, int row_size)
: Block(0, row_size), pos{pos}, _edges(this), _elems(this), _alive{true}, _mass{1}, _glued_to(this) {
}

void Vertex::eat(Vertex& other) {
  if (&other == this) return;
  other._alive = false;
  pos = (_mass*pos + other._mass*other.pos)/(_mass + other._mass);
  _mass += other._mass;
  for (int i = other._edges.partners().size() - 1; i >= 0; --i) pair(other._edges.partners()[i]);
  for (int i = other._elems.partners().size() - 1; i >= 0; --i) pair(other._elems.partners()[i]);
  other._mass = 0;
}

void Vertex::glue(Element_shape& to, std::vector<double> coords) {
  HEXED_ASSERT(std::size_t(to.n_dim()) == coords.size(), "wrong number of glued coordinates");
  _glued_to.pair(to._glued_verts);
  _glued_coords = coords;
}

std::vector<int> interior_dims(int n_dim, int row_size) {
  std::vector<int> dims(n_dim, row_size - 2);
  dims.push_back(3);
  return dims;
}

Boundary_interior::Boundary_interior(int n_dim, const Basis& b)
: Block(n_dim, b.row_size), _basis{&b}, _interior(interior_dims(n_dim, b.row_size)) {
}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, const Basis& b)
: Boundary_interior(1, b), _verts{this, this} {
  vertex0.pair(_verts[0]);
  vertex1.pair(_verts[1]);
  reset();
}

Mat<3> Edge::_point(std::vector<int> coords) const {
  int coord = coords[0];
  if (glued()) {
    if (_half == no) return _glued_to->_point(coords);
    else {
      Mat<3, dyn> pts(3, row_size());
      for (int c = 0; c < row_size(); ++c) pts(all, c) = _glued_to->point({c});
      return pts*basis().prolong(_half)(coord, all).transpose();
    }
  }
  if (coord ==       0) return _verts[0]->point({});
  if (coord == row_size() - 1) return _verts[1]->point({});
  return _interior(coord - 1).vector();
}

void Edge::reset() {
  for (int i = 1; i < row_size() - 1; ++i) {
    double n = basis().node(i);
    _interior(i - 1) = Array<double>(Mat<3>((1 - n)*_verts[0]->point({}) + n*_verts[1]->point({})));
  }
}

const int Edge::no = -1;

void Edge::glue(Edge& other, int half) {
  _glued_to.set(&other);
  _half = half;
}

Mat<3> Face::_point(std::vector<int> coords) const {
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    if (coords[i_dim] ==       0) return _edges[2*i_dim    ].point({coords[!i_dim]});
    if (coords[i_dim] == row_size() - 1) return _edges[2*i_dim + 1].point({coords[!i_dim]});
  }
  return _interior(coords[0] - 1)(coords[1] - 1).vector();
}

Face::Face(std::array<Vertex*, 4> verts, const Basis& b)
: Boundary_interior(2, b) {
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _edges.emplace_back(*verts[(2 - i_dim)*sign], *verts[(2 - i_dim)*sign + 1 + i_dim], basis());
    }
  }
  reset();
}

void Face::reset() {
  int rs = row_size();
  int int_sz = (rs - 2)*(rs - 2);
  int tot_sz = rs*rs;
  Mat<dyn, dyn> dmsq = basis().diff_mat()*basis().diff_mat();
  Mat_rm<> lhs_mat = Mat_rm<>::Zero(tot_sz, int_sz);
  Mat_rm<> rhs_mat = Mat_rm<>::Zero(tot_sz, 3);
  Array<double> lhs({rs, rs, rs - 2, rs - 2}, lhs_mat.data());
  Array<double> rhs({rs, rs, 3}, rhs_mat.data());
  for (int i_row = 0; i_row < rs; ++i_row) {
    for (int j_row = 1; j_row < rs - 1; ++j_row) {
      for (int col = 1; col < rs - 1; ++col) {
        lhs(i_row)(j_row)(col - 1)[j_row - 1] += dmsq(i_row, col);
        lhs(j_row)(i_row)(j_row - 1)[col - 1] += dmsq(i_row, col);
      }
      for (int col : {0, rs - 1}) {
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(    bool(col)).point({j_row});
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(2 + bool(col)).point({j_row});
      }
    }
    for (int j_row : {0, rs - 1}) {
      for (int col = 0; col < rs; ++col) {
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(    bool(j_row)).point({col});
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(2 + bool(j_row)).point({col});
      }
    }
  }
  Mat_rm<> soln = lhs_mat.fullPivHouseholderQr().solve(rhs_mat);
  _interior = soln.data();
}

int vstride(int n_dim, int i_dim) {return math::pow(2, n_dim - 1 - i_dim);}

Mat<3> Element_shape::_vertex_point(std::vector<int> coords) const {
  Mat<3> point = Mat<3>::Zero();
  for (int i_vert = 0; i_vert < math::pow(2, n_dim()); ++i_vert) {
    double weight = 1;
    for (int i_dim = 0; i_dim < n_dim(); ++i_dim) {
      bool sign = i_vert/vstride(n_dim(), i_dim)%2;
      weight *= !sign + math::sign(sign)*_basis->node(coords[i_dim]);
    }
    point += weight*_verts[i_vert]->pos;
  }
  return point;
}

Mat<3> Element_shape::_point(std::vector<int> coords) const {
  Mat<3> point = _vertex_point(coords);
  if (_i_bf != Mesh_blocks::no_face) {
    int sign = _i_bf%2;
    int i_dim = _i_bf/2;
    double interp_coef = !sign + math::sign(sign)*_basis->node(coords[i_dim]);
    coords[i_dim] = sign*(row_size() - 1);
    Mat<3> uncorrected = _vertex_point(coords);
    coords.erase(coords.begin() + i_dim);
    point += interp_coef*(_bf->point(coords) - uncorrected);
  }
  return point;
}

Element_shape::Element_shape(int nd, const Basis& b)
: Block(nd, b.row_size), _basis{&b}, _i_bf{6}, _glued_verts(this) {
  for (int i_vert = 0; i_vert < math::pow(2, nd); ++i_vert) _verts.emplace_back(this);
}

int i_edge(Connection_direction dir, int side, int i_bf) {
  return 2*(dir.i_dim[side] > 3 - dir.i_dim[side] - i_bf/2) + dir.face_sign[side];
}

#define ASSERT_CON_DIMS(dir, other) \
  HEXED_ASSERT((other)._bf, "connection expected subordinate element to have boundary face"); \
  if (dir.i_dim[0] == dir.i_dim[1]) { \
    HEXED_ASSERT(_i_bf == (other)._i_bf, "boundary face mismatch on same-dim connection"); \
  } else { \
    HEXED_ASSERT(_i_bf == 2*dir.i_dim[1] + dir.face_sign[1], "boundary face mismatch on coarse element"); \
    HEXED_ASSERT((other)._i_bf == 2*dir.i_dim[0] + dir.face_sign[0], "boundary face mismatch in on fine element"); \
  } \

void Element_shape::connect(Element_shape& other, Connection_direction dir) {
  HEXED_ASSERT(other.n_dim() == n_dim(), "attempt to connect elements with different dimensionality");
  HEXED_ASSERT(other._basis == _basis, "attempt to connect elements with different basis");
  auto inds = vertex_inds(n_dim(), dir);
  for (int i_vert = 0; i_vert < math::pow(2, n_dim() - 1); ++i_vert) {
    vertex(inds[0][i_vert]).eat(other.vertex(inds[1][i_vert]));
  }
  if (n_dim() == 3 && _bf && _i_bf/2 != dir.i_dim[0]) {
    ASSERT_CON_DIMS(dir, other);
    other._sf->edge(i_edge(dir, 1, other._i_bf)).glue(_sf->edge(i_edge(dir, 0, _i_bf)));
  }
}

void Element_shape::connect(std::vector<Element_shape*> others, Connection_direction dir) {
  HEXED_ASSERT(others.size() == math::pow(std::size_t(2), n_dim() - 1), "wrong number of fine elements");
  for (Element_shape* other : others) {
    HEXED_ASSERT(other, "fine element pointer is null");
    HEXED_ASSERT(other->n_dim() == n_dim(), "attempt to connect elements with different dimensionality");
    HEXED_ASSERT(other->_basis == _basis, "attempt to connect elements with different basis");
  }
  auto inds = vertex_inds(n_dim(), dir);
  auto face_inds = face_vertex_inds(n_dim(), dir);
  int nv = math::pow(2, n_dim() - 1);
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    vertex(inds[0][i_vert]).eat(others[face_inds[i_vert]]->vertex(inds[1][i_vert]));
    for (int j_vert = 0; j_vert < nv; ++j_vert) if (j_vert != i_vert) {
      std::vector<double> coords(n_dim());
      bool do_it = true;
      for (int i_dim = 0, face_dim = 0; i_dim < n_dim(); ++i_dim) {
        if (i_dim == dir.i_dim[0]) coords[i_dim] = dir.face_sign[0];
        else {
          int vs = vstride(n_dim() - 1, face_dim++);
          int c = i_vert/vs%2 + j_vert/vs%2;
          do_it = do_it && !(c == 1 &&    others[face_inds[i_vert - i_vert/vs%2*vs]]
                                       == others[face_inds[i_vert + (1 - i_vert/vs%2)*vs]]);
          coords[i_dim] = .5*c;
        }
      }
      if (do_it) others[face_inds[i_vert]]->vertex(inds[1][j_vert]).glue(*this, coords);
    }
  }
  if (n_dim() == 3 && _bf && _i_bf/2 != dir.i_dim[0]) {
    Edge& edge = _sf->edge(i_edge(dir, 0, _i_bf));
    int edge_dim = _i_bf/2 > 3 - dir.i_dim[0] - _i_bf/2;
    int strides [2] {vstride(2, edge_dim), vstride(2, !edge_dim)};
    Array<Element_shape*> to_glue({2}, [&](int i){return others[face_inds[_i_bf%2*strides[0] + i*strides[1]]];});
    auto glue = [&](int i_glue, int i_half) {
      ASSERT_CON_DIMS(dir, *to_glue[i_glue]);
      to_glue[i_glue]->_sf->edge(i_edge(dir, 1, to_glue[i_glue]->_i_bf)).glue(edge, i_half);
    };
    if (to_glue[0] == to_glue[1]) glue(0, Edge::no);
    else for (int i_glue = 0; i_glue < 2; ++i_glue) glue(i_glue, i_glue);
  }
}

const int Mesh_blocks::no_face = -1;

Mesh_blocks::Mesh_blocks(int nd, const Basis& b)
: n_dim{nd}, basis{b} {
}

template <typename T>
Sequence<T&> purge_fetch(std::vector<T>& vec) {
  std::erase_if(vec, [](T& t){return !t.alive();});
  return Sequence<T&>::vector_view(vec);
}

Sequence<Vertex&> Mesh_blocks::interior_verts() {return purge_fetch(_interior_verts);}
Sequence<Vertex&> Mesh_blocks::boundary_verts() {return purge_fetch(_boundary_verts);}
Sequence<Edge&> Mesh_blocks::edges_2d() {return purge_fetch(_edges_2d);}
Sequence<Face&> Mesh_blocks::faces_3d() {return purge_fetch(_faces_3d);}

std::unique_ptr<Element_shape> Mesh_blocks::create_element(Mat<3> pos, double size, int boundary_face) {
  std::unique_ptr<Element_shape> ptr(new Element_shape(n_dim, basis));
  ptr->_i_bf = boundary_face;
  int nv = math::pow(2, n_dim);
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    auto vec = &_interior_verts;
    if (boundary_face != no_face) {
      if ((i_vert/vstride(n_dim, boundary_face/2))%2 == boundary_face%2) vec = &_boundary_verts;
    }
    Mat<3> p = pos;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) p(i_dim) += i_vert/vstride(n_dim, i_dim)%2*size;
    vec->emplace_back(p, basis.row_size);
    vec->back().pair(ptr->_verts[i_vert]);
  }
  if (boundary_face != no_face) {
    int i_dim = boundary_face/2;
    int sign = boundary_face%2;
    if (n_dim == 2) {
      int vert0 = sign*vstride(2, i_dim);
      _edges_2d.emplace_back(ptr->vertex(vert0), ptr->vertex(vert0 + vstride(2, !i_dim)), basis);
      ptr->_bf.set(&_edges_2d.back());
    } else if (n_dim == 3) {
      std::array<Vertex*, 4> verts;
      for (int i_vert = 0; i_vert < 4; ++i_vert) verts[i_vert] = &_boundary_verts.end()[i_vert - 4];
      _faces_3d.emplace_back(verts, basis);
      ptr->_bf.set(&_faces_3d.back());
      ptr->_sf.set(&_faces_3d.back());
    }
  }
  return ptr;
}

}
