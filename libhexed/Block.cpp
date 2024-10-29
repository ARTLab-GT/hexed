#include <hexed/Block.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/vertex_inds.hpp>
#include <hexed/Mesh_assessment.hpp>

namespace hexed::next {

int vstride(int n_dim, int i_dim) {return math::pow(2, n_dim - 1 - i_dim);}

Array<double> Block::points() const {
  // construct an Array with the correct shape
  std::vector<Int> shape {3};
  for (int i_dim = 0; i_dim < _n_dim; ++i_dim) shape.push_back(_row_size);
  Array<double> pts(shape);
  // populate array
  std::vector<int> inds(_n_dim);
  for (int i_point = 0; i_point < pts.stride(0); ++i_point) {
    for (int i_dim = 0; i_dim < _n_dim; ++i_dim) inds[i_dim] = (i_point/pts.stride(1 + i_dim))%_row_size;
    auto pt = point(inds);
    for (int i_dim = 0; i_dim < 3; ++i_dim) pts(i_dim)[i_point] = pt(i_dim);
  }
  return pts;
}

void Block::visualize(std::string format, std::string file_name, double time) const {
  Sequence<const Block&> seq([this](std::size_t)->const Block& {return *this;}, []()->std::size_t{return 1;});
  visualize(format, file_name, seq, time);
}

void Block::visualize(std::string format, std::string file_name, next::Sequence<const Block&> blocks, double time) {
  // construct the Visualizer
  int block_dim = blocks.empty() ? 1 : blocks[0]._n_dim;
  auto visualizer = Visualizer::create(format, 3, block_dim, file_name, {}, time, Visualizer::block);
  // write each block via the Visualizer
  for (const Block& block : blocks) {
    visualizer->write_block(block.points(), Array<double>({}));
  }
}

Mat<3> Block::point(const std::vector<int>& node_coords) const {
  #ifdef DEBUG
  HEXED_ASSERT(int(node_coords.size()) == _n_dim, "wrong number of node coordinates");
  for (int coord : node_coords) {
    HEXED_ASSERT(0 <= coord && coord < _row_size, format_str(1000, "node index %i is out of bounds", coord));
  }
  #endif
  return _point(node_coords);
}

Mat<3> Block::point(int i_point) const {
  #ifdef DEBUG
  HEXED_ASSERT(0 <= i_point && i_point < math::pow(_row_size, _n_dim),
               format_str(1000, "node index %i is out of bounds", i_point));
  #endif
  std::vector<int> node_coords(_n_dim);
  for (int i_dim = _n_dim - 1, stride = 1; i_dim >= 0; --i_dim, stride *= _row_size) {
    node_coords[i_dim] = (i_point/stride)%_row_size;
  }
  return point(node_coords);
}

Vertex::Vertex(Mat<3> pos, int row_size)
: Block(0, row_size)
, snapped_edge{-1}
, snapped_endpoint{-1}
, _pos{pos}
, _update{Mat<3>::Zero()}
, _target{Mat<3>::Zero()}
, _has_target{false}
, _edges(this)
, _elems(this)
, _glued_to(this)
, _shadowed(this)
, _shadows(this)
, _shared_value{0}
{}

Vertex::~Vertex() {
  for (auto v : _shadows.theirs()) v->_pos = point({});
}

double Vertex::nominal_size() const {
  double nom_sz = 0;
  for (auto elem : _elems.theirs()) {
    HEXED_ASSERT("elem", "element is null");
    nom_sz = std::max(nom_sz, elem->nominal_size());
  }
  return nom_sz;
}

void Vertex::shadow(Vertex& that) {
  that._pos = _pos = .5*(that.point({}) + point({}));
  HEXED_ASSERT(that._shadowed.get() != this, "two `Vertex`s cannot shadow each other");
  HEXED_ASSERT(!_shadowed || !that._shadowed, "one of the vertices must not already be shadowing");
  if (_shadowed) that._shadowed.pair(_shadows);
  else _shadowed.pair(that._shadows);
}

void Vertex::eat(Vertex& that) {
  HEXED_ASSERT(alive() && that.alive(), "both vertices must be alive (at least at the start...)");
  if (&that == this) return;
  // compute averaged position
  Int sz [2] {_elems.partners().size(), that._elems.partners().size()};
  _pos = (sz[0]*point({}) + sz[1]*that.point({}))/(sz[0] + sz[1]);
  // steal pointers
  for (Int i = that._edges.partners().size() - 1; i >= 0; --i) pair(that._edges.partners()[i]);
  for (Int i = that._elems.partners().size() - 1; i >= 0; --i) pair(that._elems.partners()[i]);
  record.insert(record.end(), that.record.begin(), that.record.end());
  if (that.snapped_endpoint >= 0 && snapped_endpoint < 0) {
    snapped_edge = that.snapped_edge;
    snapped_endpoint = that.snapped_endpoint;
  }
}

void Vertex::glue(Element_shape& to, std::vector<double> coords) {
  HEXED_ASSERT(std::size_t(to.n_dim()) == coords.size(), "wrong number of glued coordinates");
  _glued_to.pair(to._glued_verts);
  _glued_coords = coords;
}

void Vertex::calc_relax() {
  _update = _desired_pos() - _pos;
}

void Vertex::apply_relax() {
  if (_shadowed || glued()) return;
  Mat<3> u = _update;
  for (auto s : _shadows.theirs()) u += s->_update;
  _pos += u/(1 + _shadows.theirs().size());
}

double Vertex::badness(Mat<3> proposed_pos) const {
  Mat<3> des_pos = _desired_pos();
  for (auto s : _shadows.theirs()) des_pos += s->_desired_pos();
  des_pos /= 1 + _shadows.theirs().size();
  return (proposed_pos - des_pos).norm();
}

void Vertex::set_pos(Mat<3> p) {
  if (!glued()) {
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      #pragma omp atomic write
      _pos(i_dim) = p(i_dim);
    }
  }
}

void Vertex::reset_pos() {
  Mat<3> p = Mat<3>::Zero();
  int n = 0;
  for (auto elem : _elems.theirs()) if (elem) if (!elem->glued()) {
    p += elem->nominal_position(_get_index(*elem));
    ++n;
  }
  if (n) _pos = p/n;
}

constexpr double jacobian_tolerance = 1e-2;

double compute_badness(double value, double target, double lower_bound) {
  return math::pow((value - target)/(value - lower_bound*target), 2);
}

double deriv_badness(double value, double target, double lower_bound) {
  double num = value - target;
  double denom = value - lower_bound*target;
  return 2*num/denom*(1/denom - num/denom/denom);
}

bool Vertex::mobile() const {
  bool m = false;
  for (auto elem : _elems.theirs()) if (elem) m = m || (!elem->glued() && elem->deformed);
  m = m && !glued();
  return m;
}

Vertex::_Optimization_state Vertex::_compute_state() {
  set_pos(point({}));
  _Optimization_state state;
  state.feasible = true;
  state.objective = 0;
  state.gradient.setZero();
  int nd = _elems.theirs()[0]->n_dim();
  int nv = math::pow(2, nd);
  for (auto elem : _elems.theirs()) {
    HEXED_ASSERT(elem, "element is null");
    if (elem->glued()) continue;
    int i_this = _get_index(*elem);
    double ns = elem->nominal_size();
    Mat<3, dyn> verts(3, nv);
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      verts(all, i_vert) = elem->vertex(i_vert).point({});
    }
    Sequence<Mat<3>> vert_seq {
      [&](Int i_vert)->Mat<3> {return verts(all, i_vert);},
      [&]()->Int {return nv;},
    };
    std::vector<int> i_those(nd + 1);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      i_those[i_dim] = i_this - math::sign(i_this/vstride(nd, i_dim)%2)*vstride(nd, i_dim);
      HEXED_ASSERT(i_those[i_dim] >= 0 && i_those[i_dim] < nv, "`i_those` out of bounds");
    }
    i_those[nd] = i_this;
    for (int i_that : i_those) {
      Mesh_assessment ma(vert_seq, i_that, i_this);
      state.feasible = state.feasible && ma.orthogonality > jacobian_tolerance;
      for (int i_dim = 0; i_dim < nd; ++i_dim) state.feasible = state.feasible && ma.edge_lengths(i_dim) > jacobian_tolerance*ns;
      if (state.feasible) {
        state.objective += 10*compute_badness(ma.orthogonality, 1., jacobian_tolerance);
        state.gradient += 10*deriv_badness(ma.orthogonality, 1., jacobian_tolerance)*ma.grad_orth;
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          state.objective += compute_badness(ma.edge_lengths(i_dim), ns, jacobian_tolerance);
          state.gradient += deriv_badness(ma.edge_lengths(i_dim), ns, jacobian_tolerance)*ma.grad_lengths(i_dim, all).transpose();
        }
        if (!glued() && elem->vertex(i_that).glued()) {
          auto that_state = elem->vertex(i_that)._compute_state();
          state.objective += that_state.objective;
          state.gradient += .5*that_state.gradient;
          state.feasible = state.feasible && that_state.feasible;
        }
      }
    }
  }
  return state;
}

void Vertex::improve_quality() {
  auto state = _compute_state();
  double ns = nominal_size();
  HEXED_ASSERT(state.feasible, "Vertex state violates quality criteria.");
  //if (state.gradient.norm()*ns > 1e-4*state.objective) {
  if (true) {
    _Optimization_state new_state;
    double step_sz = .1*ns;
    Mat<3> step_dir = -state.gradient.normalized();
    Mat<3> orig_pos = _point({});
    do {
      if (step_sz < 1e-12*ns) {
        set_pos(orig_pos);
        std::cout << "  step rejected in `improve_quality`" << std::endl;
        break;
      }
      set_pos(orig_pos + step_sz*step_dir);
      new_state = _compute_state();
      step_sz /= 2;
    } while (!(new_state.feasible && new_state.objective < state.objective));
  }
}

void Vertex::move_toward(std::function<Mat<3>(Mat<3>)> get_target) {
  Mat<3> orig_pos = _point({});
  auto state = _compute_state();
  double ns = nominal_size();
  HEXED_ASSERT(state.feasible, "Vertex state violates quality criteria.");
  Mat<3> improve_dir = -state.gradient.normalized();
  double improve_sz = .1*ns;
  Mat<3> target = get_target(orig_pos);
  Mat<3> snap_dir = target - orig_pos;
  double snap_sz = snap_dir.norm();
  double orig_dist = snap_sz;
  snap_dir /= snap_sz;
  _Optimization_state new_state;
  do {
    if (improve_sz < 1e-12*ns) {
      _pos = orig_pos;
      std::cout << "  improvement step rejected in `move_toward`" << std::endl;
      break;
    }
    Mat<3> improved_pos = orig_pos + improve_sz*improve_dir;
    _pos = improved_pos;
    new_state = _compute_state();
    improve_sz /= 2;
    target = get_target(improved_pos);
    if (new_state.feasible && new_state.objective < state.objective) {
      do {
        if (snap_sz < 1e-14*ns) {
          _pos = improved_pos;
          std::cout << "    snapping step rejected in `move_toward`" << std::endl;
          break;
        }
        _pos = improved_pos + snap_sz*snap_dir;
        new_state = _compute_state();
        snap_sz /= 2;
      } while (!(new_state.feasible && new_state.objective < 10*state.objective));
    }
  } while (!(new_state.feasible && (target - _pos).norm() <= orig_dist + 1e-8*ns));
}

void Vertex::set_target(Mat<3> p) {
  _target = p;
  _has_target = true;
}

void Vertex::set_target() {
  _has_target = false;
}

double Vertex::quality() {
  auto state = _compute_state();
  HEXED_ASSERT(state.feasible, "Vertex state violates quality criteria.");
  return state.objective;
}

int Vertex::n_elements() const {
  int n = 0;
  for (auto elem : _elems.theirs()) n += bool(elem);
  return n;
}

#define NEIGHBORS(CONST) \
std::vector<CONST Vertex*> Vertex::neighbors() CONST { \
  std::vector<CONST Vertex*> n; \
  for (auto elem : _elems.theirs()) { \
    HEXED_ASSERT(_elems.theirs()[0], "element is null"); \
    int nv = math::pow(2, elem->n_dim()); \
    int i_this = _get_index(*elem); \
    for (int i_vert = 0; i_vert < nv; ++i_vert) { \
      CONST Vertex* vert = &elem->vertex(i_vert); \
      for (int stride = 1; stride < nv; stride *= 2) { \
        if (i_this - i_vert == stride*(2*(i_this/stride%2) - 1) && std::find(n.begin(), n.end(), vert) == n.end()) { \
          n.push_back(vert); \
        } \
      } \
    } \
  } \
  return n; \
}
NEIGHBORS()
NEIGHBORS(const)
#undef NEIGHBORS

Vertex::Shared_value::Shared_value(Vertex& vert) : _vert{vert} {
  if (!_vert.glued()) _acquire.emplace(_vert._shared_value_lock);
}

double Vertex::Shared_value::get() const {
  if (_acquire) return _vert._shared_value;
  HEXED_ASSERT(_vert.glued(), "The glued status of the vertex changed since constructing the `Shared_value`.")
  double value = 0;
  for (int i_vert = 0; i_vert < math::pow(2, _vert._glued_to->n_dim()); ++i_vert) {
    double interp = 1.;
    bool skip = false;
    for (int i_dim = 0; i_dim < _vert._glued_to->n_dim(); ++i_dim) {
      int sign = i_vert/vstride(_vert._glued_to->n_dim(), i_dim)%2;
      skip = skip || (_vert._glued_coords[i_dim] == !sign);
      interp *= !sign + math::sign(sign)*_vert._glued_coords[i_dim];
    }
    if (!skip) value += interp*Shared_value(_vert._glued_to->vertex(i_vert)).get();
  }
  return value;
}

void Vertex::Shared_value::set(double value) {
  if (_acquire) _vert._shared_value = value;
}

Mat<3> Vertex::_point(const std::vector<int>&) const {
  // usually, the vertex will not be glued or a shadow and we can just return the `_pos`
  if (_shadowed) return _shadowed.value().point({});
  if (!_glued_to) {
    Mat<3> p;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      #pragma omp atomic read
      p(i_dim) = _pos(i_dim);
    }
    return p;
  }
  return _glued_to.value().interpolate(_glued_coords);
}

Mat<3> Vertex::_desired_pos() const {
  Mat<3> des_pos = Mat<3>::Zero();
  HEXED_ASSERT(alive(), "`Vertex` must be `alive()` to compute optimize postion");
  HEXED_ASSERT(_elems.theirs()[0], "element is null");
  #if 1
  auto n = neighbors();
  for (const Vertex* vert : n) {
    des_pos += vert->point({});
  }
  des_pos = .1*_pos + .9*des_pos/n.size();
  #else
  int nd = _elems.theirs()[0]->n_dim();
  int nv = math::pow(2, nd);
  double tot_sz = 0;
  for (auto elem : _elems.theirs()) {
    HEXED_ASSERT(elem, "element is null");
    if (elem->glued()) continue;
    double nom_sz = elem->nominal_size();
    int i_this = _get_index(*elem);
    Mat<3, dyn> verts(3, nv);
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      // we can use `_pos` because `Mesh_blocks` just set that to `point({})`
      verts(all, i_vert) = elem->vertex(i_vert)._pos;
    }
    HEXED_ASSERT(i_this >= 0, "`this` does not appear to be a vertex of `elem`!");
    if (!elem->deformed) return elem->nominal_position(i_this);
    std::vector<int> coords(nd);
    for (int i_dim = 0; i_dim < nd; ++i_dim) coords[i_dim] = i_this/vstride(nd, i_dim)%2*vstride(nd, i_dim);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      Mat<3, 2> edges;
      int opposite = i_this + (vstride(nd, i_dim) - 2*coords[i_dim]);
      bool degenerate = false;
      for (int i_edge = 0; i_edge < 2; ++i_edge) {
        if (i_edge < nd - 1) {
          int j_dim = (i_dim + i_edge + 1)%nd;
          int start = opposite - coords[j_dim];
          edges(all, i_edge) = verts(all, start + vstride(nd, j_dim)) - verts(all, start);
          if (edges(all, i_edge).norm() < 1e-2*nom_sz) degenerate = true;
        } else edges(all, i_edge) = math::sign(!i_dim)*nom_sz*Mat<3>::Unit(2);
      }
      if (!degenerate) {
        tot_sz += 1/nom_sz;
        des_pos += (verts(all, opposite) + math::sign(coords[i_dim])/nom_sz*edges(all, 0).cross(edges(all, 1)))
                   /nom_sz;
      }
    }
  }
  if (tot_sz == 0) return _pos;
  des_pos = .9*des_pos/tot_sz + .1*_pos;
  #endif
  return des_pos;
}

int Vertex::_get_index(const Element_shape& elem) const {
  int i_this = -1;
  int nv = math::pow(2, elem.n_dim());
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    if (&elem.vertex(i_vert) == this) i_this = i_vert;
  }
  HEXED_ASSERT(i_this >= 0, "`this` is not a vertex of `elem`");
  return i_this;
}

std::vector<Int> interior_dims(int n_dim, int row_size) {
  std::vector<Int> dims(n_dim, row_size - 2);
  dims.push_back(3);
  return dims;
}

Boundary_block::Boundary_block(int n_dim, const Basis& b)
: Block(n_dim, b.row_size)
, _interior(interior_dims(n_dim, b.row_size))
, _basis{&b}
, _elem(this)
{}

Mat<3> Edge::_point(const std::vector<int>& coords) const {
  int coord = coords[0];
  if (glued()) {
    if (_glued_reverse) {
      coord = row_size() - 1 - coord;
    }
    if (_half == no) {
      return _glued_to.value()._point({coord});
    } else {
      Mat<3, dyn> pts(3, row_size());
      for (int c = 0; c < row_size(); ++c) pts(all, c) = _glued_to.value().point({c});
      return pts*basis().prolong(_half)(coord, all).transpose();
    }
  }
  if (coord ==       0) return _verts[0].value().point({});
  if (coord == row_size() - 1) return _verts[1].value().point({});
  return _interior(coord - 1).vector();
}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, const Basis& b)
: Boundary_block(1, b)
, _verts{this, this}
, _glued_to(this)
, _glued(this)
, _half{no}
, _glued_reverse{false}
{
  vertex0.pair(_verts[0]);
  vertex1.pair(_verts[1]);
  reset();
}

void Edge::reset() {
  for (int i = 1; i < row_size() - 1; ++i) {
    double n = basis().node(i);
    _interior(i - 1) = Array<double>(Mat<3>((1 - n)*_verts[0].value().point({}) + n*_verts[1].value().point({})));
  }
}

const int Edge::no = -1;

void Edge::glue(Edge& other, int half, bool reverse) {
  HEXED_ASSERT(&other != this, "cannot glue an edge to itself");
  HEXED_ASSERT(!other._glued_to, "Cascading edge gluing is forbidden (in order to catch algorithmic bugs).");
  _glued_to.pair(other._glued);
  _half = half;
  _glued_reverse = reverse;
}

bool Edge::glued() const {
  if (_glued_to) return _glued_to->alive();
  else return false;
}

std::vector<Element_shape*> Edge::contacted_elements() {
  std::vector<Element_shape*> elems;
  if (element()) elems.push_back(element());
  if (glued()) elems.push_back(_glued_to->element());
  for (Edge* edge : _glued.theirs()) if (edge) if (edge->alive()) elems.push_back(edge->element());
  return elems;
}

Mat<3> Face::_point(const std::vector<int>& coords) const {
  // if the point is on the boundary of the node array, forward to one of the edges
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    if (coords[i_dim] ==              0) return _edges[2*i_dim    ].point({coords[!i_dim]});
    if (coords[i_dim] == row_size() - 1) return _edges[2*i_dim + 1].point({coords[!i_dim]});
  }
  // otherwise, return an interior node
  return _interior(coords[0] - 1)(coords[1] - 1).vector();
}

Face::Face(std::array<Vertex*, 4> verts, const Basis& b) : Boundary_block(2, b) {
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      _edges.emplace_back(*verts[(2 - i_dim)*sign], *verts[(2 - i_dim)*sign + 1 + i_dim], basis());
    }
  }
  reset();
}

void Face::reset() {
  int rs = row_size();
  int interior_sz = (rs - 2)*(rs - 2);
  int total_sz = rs*rs;
  // 1D second derivative matrix
  Mat<dyn, dyn> dmsq = basis().diff_mat()*basis().diff_mat();
  // the LHS matrix computes the influence of each interior point on the Laplacian of all points
  Mat_rm<> lhs_mat = Mat_rm<>::Zero(total_sz, interior_sz);
  // each column of the RHS matrix is the Laplacian of one physical coordinate at all nodes
  // when the interior nodes are all zero
  Mat_rm<> rhs_mat = Mat_rm<>::Zero(total_sz, 3);
  // array views of the matrices
  Array<double> lhs({rs, rs, rs - 2, rs - 2}, lhs_mat.data());
  Array<double> rhs({rs, rs, 3}, rhs_mat.data());
  // populate arrays/matrices
  for (int i_row = 0; i_row < rs; ++i_row) {
    for (int j_row = 1; j_row < rs - 1; ++j_row) {
      // add influence of interior nodes
      for (int col = 1; col < rs - 1; ++col) {
        lhs(i_row)(j_row)(col - 1)[j_row - 1] += dmsq(i_row, col);
        lhs(j_row)(i_row)(j_row - 1)[col - 1] += dmsq(i_row, col);
      }
      // add influence of edge (not corner) nodes on boundary-normal rows
      for (int col : {0, rs - 1}) {
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(    bool(col)).point({j_row});
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(2 + bool(col)).point({j_row});
      }
    }
    // add influence of edges (including corners) nodes on other nodes on the same edge
    for (int j_row : {0, rs - 1}) {
      for (int col = 0; col < rs; ++col) {
        rhs(j_row)(i_row).vector() -= dmsq(i_row, col)*edge(    bool(j_row)).point({col});
        rhs(i_row)(j_row).vector() -= dmsq(i_row, col)*edge(2 + bool(j_row)).point({col});
      }
    }
  }
  // compute least-squares solution
  Mat_rm<> soln = lhs_mat.fullPivHouseholderQr().solve(rhs_mat);
  _interior = soln.data();
}

Mat<3> Element_shape::_vertex_point(const std::vector<int>& coords) const {
  // computes a point via order-1 interpolation between the vertices,
  // without accounting for boundary-side warping
  Mat<3> point = Mat<3>::Zero();
  for (int i_vert = 0; i_vert < math::pow(2, n_dim()); ++i_vert) {
    double weight = 1;
    bool skip = false;
    for (int i_dim = 0; i_dim < n_dim(); ++i_dim) {
      bool sign = i_vert/vstride(n_dim(), i_dim)%2;
      skip = skip || (coords[i_dim] == (row_size() - 1)*!sign);
      weight *= !sign + math::sign(sign)*_basis->node(coords[i_dim]);
    }
    if (!skip) point += weight*_verts[i_vert].value().point({});
  }
  return point;
}

Mat<3> Element_shape::_point(const std::vector<int>& coords) const {
  if (glued()) {
    Int nc = coords.size();
    std::vector<double> new_coords(nc);
    for (int i_dim = 0; i_dim < nc; ++i_dim) {
      double diff = _glued_corners[1][i_dim] - _glued_corners[0][i_dim];
      new_coords[i_dim] = _glued_corners[0][i_dim] + _basis->node(coords[i_dim])*diff;
    }
    return _glued_to->interpolate(new_coords);
  }
  // first compute point by interpolating between vertices
  Mat<3> point = _vertex_point(coords);
  // then, if `this` has a side on the boundary, adjust it to account for the actual position of the boundary nodes
  if (_i_bf != Mesh_blocks::no_face) {
    std::vector<int> c(coords);
    int sign = _i_bf%2;
    int i_dim = _i_bf/2;
    double interp_coef = !sign + math::sign(sign)*_basis->node(c[i_dim]);
    c[i_dim] = sign*(row_size() - 1);
    Mat<3> uncorrected = _vertex_point(c);
    c.erase(c.begin() + i_dim);
    point += interp_coef*(_bf.value().point(c) - uncorrected);
  }
  return point;
}

void Element_shape::_glue_edges(std::vector<Element_shape*> those) {
  if (n_dim() != 3 || !_sf) return;
  for (Element_shape* that : those) if (that->_sf) {
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      auto& edge0 = that->_sf.value().edge(i_edge);
      for (int j_edge = 0; j_edge < 4; ++j_edge) {
        auto& edge1 = _sf.value().edge(j_edge);
        for (bool reverse : {0, 1}) {
          if (&edge0.vertex(0) == &edge1.vertex(reverse) && &edge0.vertex(1) == &edge1.vertex(!reverse)) {
            edge0.glue(edge1, Edge::no, reverse);
          }
        }
      }
    }
  }
}

Element_shape::Element_shape(int nd, const Basis& b)
: Block(nd, b.row_size)
, deformed{false}
, extruded_direction{Mesh_blocks::no_face}
, is_new{false}
, _basis{&b}
, _i_bf{6}
, _bf(this)
, _boundary_edges(this)
, _glued_verts(this)
{
  for (int i_vert = 0; i_vert < math::pow(2, nd); ++i_vert) _verts.emplace_back(this);
}

Mat<3> Element_shape::interpolate(std::vector<double> ref_coords) const {
  int nd = n_dim();
  int rs = row_size();
  std::vector<Int> shape(nd + 1, rs);
  shape[0] = 3;
  Array<double> points(shape);
  points = 0.;
  for (int i_point = 0; i_point < (int)points.size()/3; ++i_point) {
    bool skip = false;
    std::vector<int> coords(nd);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      coords[i_dim] = i_point/math::pow(rs, nd - 1 - i_dim)%rs;
      skip = skip || (ref_coords[i_dim] == 0 && coords[i_dim] != 0       );
      skip = skip || (ref_coords[i_dim] == 1 && coords[i_dim] != (rs - 1));
    }
    if (!skip) {
      points.reshaped({3, whatever}).column(i_point).vector() = point(coords);
    }
  }
  Mat<3> p; // this is where we will put the computed position
  // compute the interpolation matrix
  Eigen::Map<const Mat<>> sample(ref_coords.data(), ref_coords.size());
  Mat<dyn, dyn> interp = _basis->interpolate(sample);
  // apply the interpolation matrix along each dimension to compute the desired point
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    Mat<> vec = points(i_dim).vector();
    for (int j_dim = nd - 1; j_dim >= 0; --j_dim) {
      Mat<> new_vec = math::dimension_matvec(interp(j_dim, all), vec, j_dim);
      vec = new_vec;
    }
    p(i_dim) = vec(0);
  }
  return p;
}

Mat<3> Element_shape::nominal_position(int i_vert) const {
  Mat<3> pos = _nom_pos;
  for (int i_dim = 0; i_dim < n_dim(); ++i_dim) {
    int sign = i_vert/vstride(n_dim(), i_dim)%2;
    pos(i_dim) += sign*_nom_sz;
    if (extruded_direction != Mesh_blocks::no_face && extruded_direction/2 == i_dim && extruded_direction%2 == sign) {
      pos(i_dim) -= math::sign(sign)*_nom_sz;
    }
  }
  return pos;
}

Mat<3> Element_shape::nominal_center() const {
  Mat<3> c = _nom_pos + Mat<3>::Constant(.5*_nom_sz);
  if (extruded_direction != Mesh_blocks::no_face) {
    c(extruded_direction/2) -= math::sign(extruded_direction%2)*.5*_nom_sz;
  }
  return c;
};

void Element_shape::connect(Element_shape& that, Connection_direction dir) {
  HEXED_ASSERT(that.n_dim() == n_dim(), "attempt to connect elements with different dimensionality");
  HEXED_ASSERT(that._basis == _basis, "attempt to connect elements with different basis");
  if (&that == this) return;
  // eat vertices
  auto inds = vertex_inds(n_dim(), dir);
  for (int i_vert = 0; i_vert < math::pow(2, n_dim() - 1); ++i_vert) {
    vertex(inds[0][i_vert]).eat(that.vertex(inds[1][i_vert]));
  }
  _glue_edges({&that});
}

void Element_shape::connect(std::vector<Element_shape*> those, Connection_direction dir) {
  HEXED_ASSERT(those.size() == math::pow(std::size_t(2), n_dim() - 1), "wrong number of fine elements");
  for (Element_shape* that : those) {
    HEXED_ASSERT(that, "fine element pointer is null");
    HEXED_ASSERT(that != this, "Self-connections must be conforming.");
    HEXED_ASSERT(that->n_dim() == n_dim(), "attempt to connect elements with different dimensionality");
    HEXED_ASSERT(that->_basis == _basis, "attempt to connect elements with different basis");
  }
  // eat/glue vertices
  auto inds = vertex_inds(n_dim(), dir);
  auto face_inds = face_vertex_inds(n_dim(), dir);
  int nv = math::pow(2, n_dim() - 1);
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    HEXED_ASSERT(those[face_inds[i_vert]], "foo");
    // eat the non-hanging vertices of the fine elements
    vertex(inds[0][i_vert]).eat(those[face_inds[i_vert]]->vertex(inds[1][i_vert]));
    // glue the hanging vertices
    for (int j_vert = 0; j_vert < nv; ++j_vert) if (j_vert != i_vert) {
      std::vector<double> coords(n_dim());
      bool do_it = true;
      for (int i_dim = 0, face_dim = 0; i_dim < n_dim(); ++i_dim) {
        if (i_dim == dir.i_dim[0]) coords[i_dim] = dir.face_sign[0];
        else {
          int vs = vstride(n_dim() - 1, face_dim++);
          int c = i_vert/vs%2 + j_vert/vs%2;
          do_it = do_it && !(c == 1 &&    those[face_inds[i_vert - i_vert/vs%2*vs]]
                                       == those[face_inds[i_vert + (1 - i_vert/vs%2)*vs]]);
          coords[i_dim] = .5*c;
        }
      }
      if (do_it) those[face_inds[i_vert]]->vertex(inds[1][j_vert]).glue(*this, coords);
    }
  }
  _glue_edges(those);
}

void Element_shape::glue(Element_shape& that, std::array<std::vector<double>, 2> corners) {
  if (that.glued()) {
    HEXED_ASSERT(that._glued_to.get() != this, "Mutually gluing 2 elements, which would create infinite recursion.");
  }
  _glued_to.set(&that);
  _glued_corners = corners;
}

const int Mesh_blocks::no_face = -1;

Mesh_blocks::Mesh_blocks(int nd, const Basis& b): n_dim{nd}, basis{b} {}

template <typename T>
Sequence<T&> purge_fetch(std::vector<T>& vec) {
  // purge the dead entries then return the alive ones
  std::erase_if(vec, [](T& t){return !t.alive();});
  return Sequence<T&>::vector_view(vec);
}

Sequence<Vertex&> Mesh_blocks::interior_verts() {return purge_fetch(_interior_verts);}
Sequence<Vertex&> Mesh_blocks::boundary_verts() {return purge_fetch(_boundary_verts);}
Sequence<Edge&> Mesh_blocks::edges_2d() {return purge_fetch(_edges_2d);}
Sequence<Face&> Mesh_blocks::faces_3d() {return purge_fetch(_faces_3d);}

Sequence<Boundary_block&> Mesh_blocks::boundary_sides() {
  if (n_dim == 2) return edges_2d().cast<Boundary_block&>();
  else if (n_dim == 3) {
    faces_3d(); // calling this performs purge
    std::vector<Face>& faces = _faces_3d;
    return {
      [&faces](std::size_t index)->Boundary_block& {
        if (index < 4*faces.size()) return faces[index/4].edge(index%4);
        else return faces[index - 4*faces.size()];
      },
      [&faces](){return 5*faces.size();},
    };
  } else return Sequence<Boundary_block&>();
}

void Mesh_blocks::relax_vertices() {
  auto vs = verts();
  #pragma omp parallel for
  for (auto& vert : vs) vert.set_pos(vert.point({}));
  #pragma omp parallel for
  for (auto& vert : vs) vert.calc_relax();
  #pragma omp parallel for
  for (auto& vert : vs) vert.apply_relax();
}

Element_shape Mesh_blocks::create_element(Mat<3> pos, double size, int boundary_face) {
  // create the element
  Element_shape elem(n_dim, basis);
  elem._nom_sz = size;
  elem._nom_pos = pos;
  // create vertices for the element and connect the element's vertex pointers to it
  int nv = math::pow(2, n_dim);
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    auto vec = &_interior_verts;
    if (boundary_face != no_face) {
      if ((i_vert/vstride(n_dim, boundary_face/2))%2 == boundary_face%2) vec = &_boundary_verts;
    }
    vec->emplace_back(elem.nominal_position(i_vert), basis.row_size);
    vec->back().pair(elem._verts[i_vert]);
  }
  // if necessary, create a `Boundary_block` and connect the element's boundary side pointer to it
  elem._i_bf = boundary_face;
  if (boundary_face != no_face) {
    int i_dim = boundary_face/2;
    int sign = boundary_face%2;
    if (n_dim == 2) { // if 2D, the `Boundary_block` is an `Edge`
      int vert0 = sign*vstride(2, i_dim);
      _edges_2d.emplace_back(elem.vertex(vert0), elem.vertex(vert0 + vstride(2, !i_dim)), basis);
      _edges_2d.back().pair(elem._bf);
    } else if (n_dim == 3) { // if 3D, the `Boundary_block` is a `Face`
      std::array<Vertex*, 4> verts;
      for (int i_vert = 0; i_vert < 4; ++i_vert) verts[i_vert] = &_boundary_verts.end()[i_vert - 4];
      _faces_3d.emplace_back(verts, basis);
      _faces_3d.back().pair(elem._bf);
      for (int i_edge = 0; i_edge < 4; ++i_edge) _faces_3d.back().edge(i_edge).pair(elem._boundary_edges);
      elem._sf.set(&_faces_3d.back());
    }
  }
  return elem;
}

}
