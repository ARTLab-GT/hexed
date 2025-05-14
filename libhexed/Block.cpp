#include <hexed/Block.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/vertex_inds.hpp>
#include <hexed/Mesh_assessment.hpp>
#include <hexed/Printer.hpp>

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

Mat<3> Block::point(const std::vector<int>& node_coords, Int recursion_depth) const {
  #ifdef DEBUG
  HEXED_ASSERT(int(node_coords.size()) == _n_dim, "wrong number of node coordinates");
  for (int coord : node_coords) {
    HEXED_ASSERT(0 <= coord && coord < _row_size, format_str(1000, "node index %i is out of bounds", coord));
  }
  #endif
  HEXED_ASSERT(recursion_depth <= max_recursion_depth, "Max recursion depth exceeded.", assert::Overflow_error);
  return _point(node_coords, recursion_depth);
}

Mat<3> Block::point(int i_point, Int recursion_depth) const {
  #ifdef DEBUG
  HEXED_ASSERT(0 <= i_point && i_point < math::pow(_row_size, _n_dim),
               format_str(1000, "node index %i is out of bounds", i_point));
  #endif
  std::vector<int> node_coords(_n_dim);
  for (int i_dim = _n_dim - 1, stride = 1; i_dim >= 0; --i_dim, stride *= _row_size) {
    node_coords[i_dim] = (i_point/stride)%_row_size;
  }
  return point(node_coords, recursion_depth);
}

Vertex::Vertex(Mat<3> pos, int row_size)
: Block(0, row_size)
, snapped_point{-1}
, snapped_edge{-1}
, snapped_endpoint{-1}
, incompatible_snap{false}
, _pos{pos}
, _orig_pos{Mat<3>::Zero()}
, _step{Mat<3>::Zero()}
, _orig_dist{0.}
, _step_sz{0.}
, _improve_failed{false}
, _improve_done{true}
, _last_step_rejected{false}
, _edges(this)
, _elems(this)
, _glued_to(this)
, _shared_value{0}
, _last_snap_failed{false}
, _sz_constraint{huge}
{}

Vertex::~Vertex() {}

double Vertex::nominal_size() const {
  double nom_sz = _sz_constraint;
  for (auto elem : _elems.theirs()) {
    HEXED_ASSERT(elem, "element is null")
    nom_sz = std::min(nom_sz, elem->nominal_size());
  }
  return nom_sz;
}

void Vertex::eat(Vertex& that) {
  HEXED_ASSERT(alive() && that.alive(), "both vertices must be alive (at least at the start...)")
  if (&that == this) return;
  // compute averaged position
  Int sz [2] {_elems.partners().size(), that._elems.partners().size()};
  // if this line throws an Overflow_error, it exits `eat` safely and there is not harm done
  set_pos((sz[0]*point({}) + sz[1]*that.point({}))/(sz[0] + sz[1]));
  // steal pointers
  for (Int i = that._edges.partners().size() - 1; i >= 0; --i) pair(that._edges.partners()[i]);
  for (Int i = that._elems.partners().size() - 1; i >= 0; --i) pair(that._elems.partners()[i]);
  record.insert(record.end(), that.record.begin(), that.record.end());
  if ((that.snapped_endpoint >= 0 && snapped_endpoint < 0) || (that.snapped_edge >= 0 && snapped_edge < 0)) {
    snapped_edge = that.snapped_edge;
    snapped_endpoint = that.snapped_endpoint;
  }
  _sz_constraint = std::min(_sz_constraint, that._sz_constraint);
}

void Vertex::glue(Element_shape& to, std::vector<double> coords) {
  HEXED_ASSERT(std::size_t(to.n_dim()) == coords.size(), "wrong number of glued coordinates");
  _glued_to.pair(to._glued_verts);
  _glued_coords = coords;
}

void Vertex::set_pos(Mat<3> p) {
  if (!glued()) {
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      #pragma omp atomic write
      _pos(i_dim) = p(i_dim);
    }
  }
}

void Vertex::add_size_constraint(double sz) {
  Lock::Set s(_shared_value_lock);
  _sz_constraint = std::min(_sz_constraint, sz);
}

void Vertex::remove_size_constraints() {
  Lock::Set s(_shared_value_lock);
  _sz_constraint = huge;
}

Mat<3> Vertex::nominal_position() const {
  Mat<3> p = Mat<3>::Zero();
  int n = 0;
  for (auto elem : _elems.theirs()) if (elem) if (!elem->glued()) {
    p += elem->nominal_position(_get_index(*elem));
    ++n;
  }
  if (n) return p/n;
  return point({});
}

bool Vertex::mobile() const {
  bool has_unglued = false;
  bool has_cartesian = false;
  for (auto elem : _elems.theirs()) if (elem) {
    has_unglued = has_unglued || !elem->glued();
    has_cartesian = has_cartesian || !(elem->glued() || elem->deformed);
  }
  return has_unglued && !has_cartesian && !glued();
}

const double ortho_tolerance = .03;
const double edge_tolerance = .02;

Vertex::_Optimization_state Vertex::_compute_state(bool include_neighbors, bool ignore, bool ignore_neighb, double extra_tol) {
  _Optimization_state state;
  _compute_state_recursive(state, 1., include_neighbors, this, ignore ? this : nullptr, ignore_neighb, extra_tol);
  return state;
}

void Vertex::_compute_state_recursive(_Optimization_state& state, double gradient_weight, bool include_neighbors,
                                      Vertex* orig_vertex, Vertex* ignore, bool ignore_neighb, double extra_tol) {
  if (!orig_vertex) orig_vertex = this;
  int nd = _elems.theirs()[0]->n_dim();
  int nv = math::pow(2, nd);
  for (Element_shape* elem : _elems.theirs()) {
    HEXED_ASSERT(elem, "element is null");
    if (elem->glued()) continue;
    int i_this = _get_index(*elem);
    Mat<3, dyn> verts(3, nv);
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      verts(all, i_vert) = elem->vertex(i_vert)._unwarped_point(ignore, ignore_neighb);
    }
    Sequence<Mat<3>> vert_seq {
      [&](Int i_vert)->Mat<3> {return verts(all, i_vert);},
      [&]()->Int {return nv;},
    };
    std::vector<int> i_those;
    if (include_neighbors) {
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        i_those.push_back(i_this - math::sign(i_this/vstride(nd, i_dim)%2)*vstride(nd, i_dim));
        HEXED_ASSERT(i_those.back() >= 0 && i_those.back() < nv, "`i_those` out of bounds");
      }
    }
    i_those.push_back(i_this);
    for (int i_that : i_those) {
      Vertex& that_vert = elem->vertex(i_that);
      double ns = that_vert.nominal_size();
      bool skip_obj = false;
      bool skip_grad = false;
      for (auto s : state.skip) if (s.elem == elem) {
        skip_obj = skip_obj || s.i == i_that;
        skip_grad = skip_grad || (s.i == i_that && s.j == i_this);
      }
      if (!skip_grad) state.skip.emplace_back(elem, i_that, i_this);
      if (state.computing_depends) {
        bool contains = false;
        for (Vertex* v : orig_vertex->_depends_on) contains = contains || &that_vert == v;
        if (!contains) orig_vertex->_depends_on.push_back(&that_vert);
      }
      Mesh_assessment ma(vert_seq, i_that, i_this);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        state.feasible = state.feasible && ma.orthogonality(i_dim) > ortho_tolerance + extra_tol;
        state.worst_ortho = std::min(state.worst_ortho, ma.orthogonality(i_dim));
        state.feasible = state.feasible && ma.edge_lengths(i_dim) > (edge_tolerance + extra_tol)*ns;
        state.worst_edge = std::min(state.worst_edge, ma.edge_lengths(i_dim)/ns);
      }
      if (state.feasible) {
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          double orth_diff = ma.orthogonality(i_dim) - ortho_tolerance;
          state.objective += (!skip_obj)*1./orth_diff;
          state.gradient += (!skip_grad)*gradient_weight*1.*(-1/orth_diff/orth_diff)
                            *ma.grad_orth(i_dim, all).transpose();
          double num = ma.edge_lengths(i_dim) - ns;
          double denom = ma.edge_lengths(i_dim) - edge_tolerance*ns;
          state.objective += (!skip_obj)*num*num/denom;
          state.gradient += (!skip_grad)*gradient_weight*(2*num/denom - num*num/(denom*denom))
                            *ma.grad_lengths(i_dim, all).transpose();
        }
        // note: the valid values for `gradient_weight` are 1., .5, and .25
        if (gradient_weight > .3 && that_vert.glued() && i_this != i_that) {
          bool coupled = false;
          for (Vertex* v : {this, orig_vertex}) {
            for (auto e : v->_elems.theirs()) {
              coupled = coupled || (!e->glued() && e == that_vert._glued_to.get());
            }
          }
          if (coupled) {
            state.has_glued_neighbor = true;
            that_vert._compute_state_recursive(state, .5*gradient_weight, true, orig_vertex, ignore, ignore_neighb, extra_tol);
          }
        }
      }
    }
  }
}

Vertex::Improve_quality_result Vertex::improve_quality() {
  return _improve_quality([](Mat<3>){return Mat<3>::Zero();}, false, false, false);
}

Vertex::Improve_quality_result Vertex::improve_quality(std::function<Mat<3>(Mat<3>)> get_target,
                                                       bool limit_direction, bool snap) {
  return _improve_quality(get_target, true, limit_direction, snap);
}

void Vertex::compute_depends() {
  _depends_on.clear();
  _Optimization_state state;
  state.computing_depends = true;
  _compute_state_recursive(state, 1., true);
  // sort to prevent deadlocks
  std::sort(_depends_on.begin(), _depends_on.end(), std::less<Vertex*>{});
}

bool Vertex::has_problem() const {
  // note that each vertex is in its own `depends_on`
  for (const Vertex* d : _depends_on) if (d->_last_snap_failed) return true;
  return false;
}

void Vertex::init_improve(std::function<Mat<3>(Mat<3>)> get_target) {
  HEXED_ASSERT(mobile(), "Only mobile vertices can be improved.")
  _orig_pos = unwarped_point();
  auto state = _compute_state();
  bool gn = false;
  for (Vertex* n : neighbors()) gn = gn || n->glued();
  HEXED_ASSERT(state.feasible, format_str(200,
               "Vertex state violates quality criteria "
               "(ortho = %e; edge = %e; coords = (%e %e %e); glued neighbor = %i).",
               state.worst_ortho, state.worst_edge, _orig_pos(0), _orig_pos(1), _orig_pos(2),
               int(gn)))
  Mat<3> grad_fd;
  bool feasible_stencil = false;
  double fd = 1e-5*nominal_size();
  for (int i = 0; i < 4 && !feasible_stencil; ++i) {
    grad_fd.setZero();
    feasible_stencil = true;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      int fd_step = 0;
      for (int sign : {-1, 1}) {
        _pos = _orig_pos;
        _pos(i_dim) += sign*fd;
        auto fd_state = _compute_state();
        if (fd_state.feasible) {
          grad_fd(i_dim) += sign*fd_state.objective;
          ++fd_step;
        } else {
          grad_fd(i_dim) += sign*state.objective;
        }
      }
      if (fd_step) grad_fd(i_dim) /= fd*fd_step;
      else feasible_stencil = false;
    }
    fd /= 10;
  }
  _pos = _orig_pos;
  _step_sz = 1.;
  _orig_dist = (get_target(_orig_pos) - _orig_pos).norm();
  _improve_done = false;
  if (feasible_stencil && grad_fd.norm()*nominal_size() > 1e-8*state.objective) {
    _step = -nominal_size()*grad_fd.normalized();
    _orig_obj = state.objective;
    _improve_failed = false;
    _last_step_rejected = false;
  } else {
    _improve_failed = true;
    _last_step_rejected = true;
  }
}

void Vertex::compute_improve(std::function<Mat<3>(Mat<3>)> get_target) {
  HEXED_ASSERT(mobile(), "Only mobile vertices can be improved.")
  if (_improve_done || _improve_failed) return;
  _step_sz /= 3;
  _pos = _orig_pos + _step_sz*_step;
  Mat<3> target = get_target(_pos);
  Mat<3> diff = target - _pos;
  double dist = diff.norm();
  if (dist > _orig_dist) {
    diff *= (dist - _orig_dist)/dist;
    _pos += diff;
  }
  if (_step_sz < 1e-10 || (_step_sz*_step + diff).norm() < 1e-3*_step_sz*_step.norm()) {
    _pos = _orig_pos;
    _improve_failed = true;
    _last_step_rejected = true;
  }
}

void Vertex::force_continue_improve() {
  _improve_done = false;
}

bool Vertex::check_improve() {
  HEXED_ASSERT(mobile(), "Only mobile vertices can be improved.")
  auto state0 = _compute_state(true, true, true);
  auto state1 = _compute_state(true, false, false, 1e-8);
  _improve_done = state0.feasible && state1.feasible && state0.objective < _orig_obj;
  return _improve_done || _improve_failed;
}

void Vertex::init_snap(std::function<Mat<3>(Mat<3>)> get_target) {
  _orig_pos = unwarped_point();
  _step = get_target(_orig_pos) - _orig_pos;
  _step_sz = 1.;
  _improve_done = false;
  _improve_failed = false;
  _last_snap_failed = false;
}

void Vertex::compute_snap() {
  _pos = _orig_pos + _step_sz*_step;
}

Vertex::Snap_result Vertex::check_snap() {
  if (_improve_failed) return {true, true, _step.norm()};
  auto state = _compute_state(true, false, false, 1e-8);
  _last_grad = state.gradient.normalized()*0.1*nominal_size();
  _improve_done = state.feasible;
  if (!_improve_done || _step_sz*_step.norm() > .2*nominal_size()*state.worst_edge) {
    _step_sz /= 3;
    _last_snap_failed = true;
  }
  _improve_failed = _step_sz < 1e-10;
  return {_improve_done, _last_snap_failed, (1 - _step_sz)*_step.norm()};
}

Int Vertex::misses = 0;
Int Vertex::tries = 0;

Vertex::Improve_quality_result Vertex::_improve_quality(std::function<Mat<3>(Mat<3>)> get_target,
                                                        bool has_target, bool limit_direction, bool snap) {
  bool miss = false;
  for (Vertex* vert : _depends_on) {
    if (!vert->_shared_value_lock.test()) {
      miss = true;
      break;
    }
  }
  if (miss) {
    #pragma omp atomic update
    ++misses;
  }
  #pragma omp atomic update
  ++tries;
  // acquire locks on all vertices we depend on to prevent simultaneous vertex motions from violating mesh quality
  std::vector<Lock::Set> locks;
  for (Vertex* vert : _depends_on) locks.emplace_back(vert->_shared_value_lock);
  _pos = unwarped_point();
  Mat<3> orig_pos = _pos;
  auto state = _compute_state();
  double ns = nominal_size();
  bool gn = false;
  for (Vertex* n : neighbors()) gn = gn || n->glued();
  HEXED_ASSERT(state.feasible, format_str(200,
               "Vertex state violates quality criteria "
               "(ortho = %e; edge = %e; coords = (%e %e %e); glued neighbor = %i).",
               state.worst_ortho, state.worst_edge, orig_pos(0), orig_pos(1), orig_pos(2),
               int(gn)))
  Mat<3> direction = -state.gradient.normalized();
  Mat<3> target = get_target(orig_pos);
  double orig_dist = (target - orig_pos).norm();
  _Optimization_state new_state = state;
  const double min_step = 1e-8*ns;
  int n_back_improve = 0;
  if (state.gradient.norm()*ns > 1e-6*state.objective) {
    //_pos = orig_pos + 1e-8*ns*direction;
    _Optimization_state test_state = _compute_state();
    #if 0
    if (std::abs(test_state.objective - state.objective + state.gradient.norm()*1e-8*ns)
        > 1e-1*std::abs(test_state.objective - state.objective)) {
      printers::warn("badGradient" + std::to_string(state.has_glued_neighbor), true);
    }
    #endif
    double step_sz = ns;
    double repeat_factor [] {10., 2., 2.};
    double end_factor [] {10., .9, 1.};
    for (int i = 0; i < 3; ++i) if (step_sz > min_step) {
      do {
        if (step_sz < min_step) {
          _pos = orig_pos;
          new_state.objective = state.objective;
          break;
        }
        step_sz /= repeat_factor[i];
        ++n_back_improve;
        _pos = orig_pos + step_sz*direction;
        Mat<3> new_target = get_target(_pos);
        double dist = (new_target - _pos).norm();
        if (dist > orig_dist) _pos += (dist - orig_dist)/dist*(new_target - _pos);
        new_state = _compute_state();
      } while (!(new_state.feasible && new_state.objective < state.objective));
      step_sz *= end_factor[i];
    }
  }
  int snap_iters = 0;
  double target_dist = 0;
  if (has_target && snap) {
    orig_pos = _pos;
    target = get_target(_pos);
    Mat<3> step = target - _pos;
    double orig_objective = new_state.objective;
    do {
      if (step.norm() < min_step) {
        _pos = orig_pos;
        new_state.objective = orig_objective;
        break;
      }
      _pos = orig_pos + step;
      target_dist = (target - _pos).norm();
      new_state = _compute_state();
      step /= 2;
      ++snap_iters;
    } while (!new_state.feasible || new_state.objective > 2*state.objective);
  }
  new_state = _compute_state();
  HEXED_ASSERT(new_state.feasible, "something changed");
  _last_snap_failed = snap_iters > 1;
  return {new_state.objective - state.objective, _last_snap_failed, target_dist};
}

bool Vertex::snap_to(Mat<3> target) {
  Mat<3> orig_pos = unwarped_point();
  set_pos(target);
  _Optimization_state state = _compute_state();
  if (!state.feasible) set_pos(orig_pos);
  return state.feasible;
}

double Vertex::quality_objective() {
  auto state = _compute_state(false);
  std::string dists;
  std::string neighb_pos;
  for (Vertex* n : neighbors()) if (n) if (n->mobile()) {
    dists += to_string(n->dijkstra_dist) + " ";
    neighb_pos += to_string(n->point({})) + " ";
  }
  HEXED_ASSERT(state.feasible, format_str(200,
               "Vertex state violates quality criteria "
               "(ortho = %e; edge = %e; coords = (%e %e %e); offset record = %s; neighbors = %s).",
               state.worst_ortho, state.worst_edge, point({})(0), point({})(1), point({})(2),
               dists.c_str(), neighb_pos.c_str()))
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
  if (!_vert.glued()) _set.emplace(_vert._shared_value_lock);
}

double Vertex::Shared_value::get() const {
  if (_set) return _vert._shared_value; // equivalent to checking if vertex is glued
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
  if (_set) _vert._shared_value = value;
}

Mat<3> Vertex::_get_pos() const {
  Mat<3> p;
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    #pragma omp atomic read
    p(i_dim) = _pos(i_dim);
  }
  return p;
}

Mat<3> Vertex::_point(const std::vector<int>&, Int recursion_depth) const {
  // usually, the vertex will not be glued and we can just return the `_pos`
  if (!_glued_to) return _get_pos();
  return _glued_to.value().interpolate(_glued_coords, recursion_depth + 1);
}

Mat<3> Vertex::unwarped_point() const {return _unwarped_point(nullptr, false);}

Mat<3> Vertex::_unwarped_point(Vertex* ignore, bool ignore_neighb) const {
  if (!_glued_to) return (((this == ignore) != ignore_neighb) && mobile()) ? _orig_pos : _get_pos();
  int nd = _glued_to->n_dim();
  int nv = math::pow(2, nd);
  int nc = _glued_coords.size();
  Mat<8, 3> vert_pos = Mat<8, 3>::Zero(); // if nd < 3, extra entries will simply have no effect
  for (int i_vert = 0; i_vert < nv; ++i_vert) {
    bool skip = false;
    for (int i_dim = 0; i_dim < nc; ++i_dim) {
      skip = skip || std::abs(_glued_coords[i_dim] - !math::row_coordinate(nd, 2, i_dim, i_vert)) < 1e-6;
    }
    if (!skip) vert_pos(i_vert, all) = _glued_to->vertex(i_vert)._unwarped_point(ignore, ignore_neighb).transpose();
  }
  Mat<3> coords = Mat<3>::Zero();
  for (int i_dim = 0; i_dim < nc; ++i_dim) coords(3 - nc + i_dim) = _glued_coords[i_dim];
  Mat<3> p;
  for (int i_dim = 0; i_dim < 3; ++i_dim) p(i_dim) = math::interp(vert_pos(all, i_dim), coords);
  return p;
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

Mat<3> Edge::_point(const std::vector<int>& coords, Int recursion_depth) const {
  int coord = coords[0];
  if (glued()) {
    if (_glued_reverse) {
      coord = row_size() - 1 - coord;
    }
    if (_half == no) {
      return _glued_to.value()._point({coord}, recursion_depth + 1);
    } else {
      Mat<3, dyn> pts(3, row_size());
      for (int c = 0; c < row_size(); ++c) pts(all, c) = _glued_to.value().point({c}, recursion_depth + 1);
      return pts*basis().prolong(_half)(coord, all).transpose();
    }
  }
  if (coord ==       0) return _verts[0].value().point({}, recursion_depth + 1);
  if (coord == row_size() - 1) return _verts[1].value().point({}, recursion_depth + 1);
  return _interior(coord - 1).vector();
}

Edge::Edge(Vertex& vertex0, Vertex& vertex1, const Basis& b)
: Boundary_block(1, b)
, snapped_edge{-1}
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

std::vector<Element_shape*> Edge::dependent_elements() {
  std::vector<Element_shape*> depend;
  if (alive() && !glued()) depend.push_back(element());
  for (Edge* e : _glued.theirs()) if (e) if (e->alive()) depend.push_back(e->element());
  return depend;
}

std::vector<int> Edge::element_coords(std::vector<int> coords) const {
  HEXED_ASSERT(element(), "must have an `element()` to call `element_coords`");
  HEXED_ASSERT(coords.size() == 1, "wrong number of edge coordinates");
  int nd = element()->n_dim();
  if (nd == 2) {
    int i_face = element()->boundary_face();
    coords.insert(coords.begin() + i_face/2, i_face%2*(row_size() - 1));
  } else if (nd == 3) {
    const Face* face = element()->boundary_face_3d();
    HEXED_ASSERT(face, "element has no face");
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      if (&face->edge(i_edge) == this) {
        coords.insert(coords.begin() + i_edge/2, i_edge%2*(row_size() - 1));
        return face->element_coords(coords);
      }
    }
    HEXED_THROW("`this` doesn't appear to be an edge of `face`");
  } else {
    HEXED_THROW("`n_dim` for an `Edge` must be 2 or 3")
  }
  return coords;
}

void Edge::reset() {
  for (int i = 1; i < row_size() - 1; ++i) {
    double n = basis().node(i);
    _interior(i - 1) = Array<double>(Mat<3>((1 - n)*_verts[0].value().point({}) + n*_verts[1].value().point({})));
  }
}

const int Edge::no = -1;

void Edge::glue(Edge& other, int half, bool reverse) {
  HEXED_ASSERT(&other != this, "cannot glue an edge to itself")
  HEXED_ASSERT(!other._glued_to, "Cascading edge gluing is forbidden (in order to catch algorithmic bugs)."
                                 " pos: " + to_string(vertex(0).unwarped_point())
                                          + to_string(vertex(1).unwarped_point()) +
                                 " pos: " + to_string(other.vertex(0).unwarped_point())
                                          + to_string(other.vertex(1).unwarped_point()))
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

int Edge::glued_half() const {
  HEXED_ASSERT(glued(), "`glued_half` is only relevant if the edge is glued")
  return _half;
}

Mat<3> Face::_point(const std::vector<int>& coords, Int recursion_depth) const {
  // if the point is on the boundary of the node array, forward to one of the edges
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    if (coords[i_dim] ==              0) return _edges[2*i_dim    ].point({coords[!i_dim]}, recursion_depth + 1);
    if (coords[i_dim] == row_size() - 1) return _edges[2*i_dim + 1].point({coords[!i_dim]}, recursion_depth + 1);
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

std::vector<Element_shape*> Face::dependent_elements() {
  std::vector<Element_shape*> depend;
  if (alive()) depend.push_back(element());
  return depend;
}

std::vector<int> Face::element_coords(std::vector<int> coords) const {
  HEXED_ASSERT(element(), "must have an `element()` to call `element_coords`");
  HEXED_ASSERT(coords.size() == 2, "wrong number of edge coordinates");
  int i_face = element()->boundary_face();
  coords.insert(coords.begin() + i_face/2, i_face%2*(row_size() - 1));
  return coords;
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

Mat<3> Element_shape::_vertex_point(const std::vector<int>& coords, Int recursion_depth) const {
  // computes a point via order 1 interpolation between the vertices,
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
    if (!skip) point += weight*_verts[i_vert].value().point({}, recursion_depth);
  }
  return point;
}

Mat<3> Element_shape::_point(const std::vector<int>& coords, Int recursion_depth) const {
  if (glued()) {
    Int nc = coords.size();
    std::vector<double> new_coords(nc);
    for (int i_dim = 0; i_dim < nc; ++i_dim) {
      double diff = _glued_corners[1][i_dim] - _glued_corners[0][i_dim];
      new_coords[i_dim] = _glued_corners[0][i_dim] + _basis->node(coords[i_dim])*diff;
    }
    return _glued_to->interpolate(new_coords, recursion_depth + 1);
  }
  // first compute point by interpolating between vertices
  Mat<3> point = _vertex_point(coords, recursion_depth + 1);
  // then, if `this` has a side on the boundary, adjust it to account for the actual position of the boundary nodes
  if (_i_bf != Mesh_blocks::no_face) {
    std::vector<int> c(coords);
    int sign = _i_bf%2;
    int i_dim = _i_bf/2;
    double interp_coef = !sign + math::sign(sign)*_basis->node(c[i_dim]);
    c[i_dim] = sign*(row_size() - 1);
    Mat<3> uncorrected = _vertex_point(c, recursion_depth + 1);
    c.erase(c.begin() + i_dim);
    point += interp_coef*(_bf.value().point(c) - uncorrected);
  }
  return point;
}

void Element_shape::_glue_edges(std::vector<Element_shape*> those) {
  if (n_dim() != 3 || !_sf) return;
  for (int i_that = 0; i_that < Int(those.size()); ++i_that) {
    Element_shape* that = those[i_that];
    if (!that->_sf) continue;
    for (int i_edge = 0; i_edge < 4; ++i_edge) {
      auto& edge0 = that->_sf.value().edge(i_edge);
      for (int j_edge = 0; j_edge < 4; ++j_edge) {
        auto& edge1 = _sf.value().edge(j_edge);
        for (bool reverse : {0, 1}) {
          if (&edge0.vertex(reverse) == &edge1.vertex(0) && &edge0.vertex(!reverse) == &edge1.vertex(1)) {
            edge0.glue(edge1, Edge::no, reverse);
          } else {
            for (int j_that = i_that + 1; j_that < Int(those.size()); ++j_that) if (those[j_that]->_sf) {
              Edge* edges [2] {&edge0, &those[j_that]->_sf.value().edge(i_edge)};
              if (   &edges[ reverse]->vertex( reverse) == &edge1.vertex(0)
                  && &edges[!reverse]->vertex(!reverse) == &edge1.vertex(1)) {
                edges[ reverse]->glue(edge1, 0, reverse);
                edges[!reverse]->glue(edge1, 1, reverse);
              }
            }
          }
        }
      }
    }
  }
}

Element_shape::Element_shape(int nd, const Basis& b)
: Block(nd, b.row_size)
, deformed{false}
, for_matching{false}
, extruded_direction{Mesh_blocks::no_face}
, is_new{false}
, record{0}
, _basis{&b}
, _i_bf{6}
, _bf(this)
, _boundary_edges(this)
, _glued_verts(this)
{
  for (int i_vert = 0; i_vert < math::pow(2, nd); ++i_vert) _verts.emplace_back(this);
}

Mat<3> Element_shape::interpolate(std::vector<double> ref_coords, Int recursion_depth) const {
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
      points.reshaped({3, whatever}).column(i_point).vector() = point(coords, recursion_depth);
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

Mat<3> Element_shape::vertex_center() const {
  Mat<3> c = Mat<3>::Zero();
  for (auto& v : _verts) c += v->point({});
  c /= _verts.size();
  return c;
}

void Element_shape::connect(Element_shape& that, Connection_direction dir) {
  HEXED_ASSERT(that._basis == _basis, "attempt to connect elements with different basis");
  if (&that == this) return;
  std::vector<Element_shape*> these(math::pow(2, n_dim() - 1), this);
  std::vector<Element_shape*> those(math::pow(2, n_dim() - 1), &that);
  connect({these, those}, dir);
  _glue_edges({&that});
}

void Element_shape::connect(std::vector<Element_shape*> those, Connection_direction dir) {
  for (Element_shape* that : those) {
    HEXED_ASSERT(that != this, "Self-connections must be conforming.");
  }
  std::vector<Element_shape*> these(math::pow(2, n_dim() - 1), this);
  connect({these, those}, dir);
  _glue_edges(those);
}

void Element_shape::connect(std::array<std::vector<Element_shape*>, 2> elems, Connection_direction dir) {
  int nd = math::log(2, elems[0].size()) + 1;
  int nv = math::pow(2, nd - 1);
  for (int i_side = 0; i_side < 2; ++i_side) {
    HEXED_ASSERT(Int(elems[i_side].size()) == nv, "wrong number of elements supplied");
    for (Element_shape* elem : elems[i_side]) {
      HEXED_ASSERT(elem, "Element pointer is null.");
      HEXED_ASSERT(elem->n_dim() == nd, "Element dimensionality does not match number of elements supplied");
    }
  }
  std::array<std::vector<int>, 2> face_inds;
  std::array<std::vector<int>, 2> inds;
  for (int i_side = 0; i_side < 2; ++i_side) {
    int rotation_sign = math::sign(dir.face_sign[i_side] == dir.face_sign[0]);
    Connection_direction side_dir {{dir.i_dim[i_side], dir.i_dim[!i_side]},
                                   {dir.face_sign[i_side], dir.face_sign[!i_side]},
                                   dir.rotate*rotation_sign};
    face_inds[i_side] = face_vertex_inds(nd, side_dir);
    inds[i_side] = vertex_inds(nd, side_dir)[0];
  }
  for (int i_elem = 0; i_elem < nv; ++i_elem) {
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      for (int i_side = 0; i_side < 2; ++i_side) {
        bool redundant = false;
        bool glue = false;
        std::vector<double> coords(nd);
        coords[dir.i_dim[i_side]] = dir.face_sign[i_side];
        for (int i_dim = 0; i_dim < nd - 1; ++i_dim) {
          int vert_coord = math::row_coordinate(nd - 1, 2, i_dim, i_vert);
          int elem_coord = math::row_coordinate(nd - 1, 2, i_dim, i_elem);
          int diff = vert_coord - elem_coord;
          int coord_dim = i_dim + (i_dim >= dir.i_dim[i_side]);
          coords[coord_dim] = vert_coord;
          if (diff != 0) {
            if (elems[!i_side][face_inds[i_side][i_elem + diff*math::stride(nd - 1, 2, i_dim)]]
                == elems[!i_side][face_inds[i_side][i_elem]]) {
              redundant = true;
            } else if (elems[i_side][i_elem + diff*math::stride(nd - 1, 2, i_dim)] == elems[i_side][i_elem]) {
              glue = true;
              coords[coord_dim] = .5;
            }
          }
        }
        if (!redundant) {
          Vertex& vert = elems[!i_side][face_inds[i_side][i_elem]]->vertex(inds[!i_side][face_inds[i_side][i_vert]]);
          if (glue) {
            vert.glue(*elems[i_side][i_elem], coords);
          } else {
            vert.eat(elems[i_side][i_elem]->vertex(inds[i_side][i_vert]));
          }
        }
      }
    }
  }
}

void Element_shape::glue(Element_shape& that, std::array<std::vector<double>, 2> corners) {
  if (that.glued()) {
    HEXED_ASSERT(that._glued_to.get() != this, "Mutually gluing 2 elements, which would create infinite recursion.");
  }
  _glued_to.set(&that);
  _glued_corners = corners;
}

void Element_shape::destroy_boundary_face() {
  _bf.unpair();
  _boundary_edges.clear();
  _sf.set();
  _i_bf = Mesh_blocks::no_face;
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
