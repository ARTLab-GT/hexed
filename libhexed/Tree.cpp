#include <queue>
#include <hexed/Tree.hpp>
#include <hexed/Row_index.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

std::array<std::vector<Element*>, 2> Tree::Connection_neighbors::elements() {
  std::array<std::vector<Element*>, 2> elems;
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (Tree* t : trees[i_side]) elems[i_side].push_back(t->elem.get());
  }
  return elems;
}

bool Tree::Connection_neighbors::valid() {
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (Tree* t : trees[i_side]) {
      if (!t) return false;
      if (!t->elem) return false;
    }
  }
  return true;
}

Tree::Tree(int nd, double root_size, Mat<> origin)
: n_dim{nd}
, elem(this)
, _root_sz{root_size}
, _ref_level{Array<int>::make_uniform({nd}, 0)}, _coords{Array<Int>::make_uniform({nd}, 0)}
, _par{nullptr}
, _children_storage()
, _face_connections(2*n_dim, nullptr)
, _fake_parents(2*n_dim, nullptr)
, _status{unprocessed}
, _is_graft{false}
{
  HEXED_ASSERT(origin.size() >= n_dim, "`origin` is too small");
  _orig = origin(Eigen::seqN(0, n_dim));
}

Tree::~Tree() {
  delete_grafts();
  for (auto c : _face_connections) HEXED_ASSERT(!c, "Attempting to destroy a tree that is still connected.")
}

Mat<> Tree::origin() const {return _orig;}
int Tree::refinement_level() const {return _ref_level.extreme(0);}

Array<int> Tree::desired_refinement_level() const {
  Array<int> rl = _ref_level.copy();
  if (elem) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) rl[i_dim] += elem->desired_refinement(i_dim);
  }
  return rl;
}

Array<int> Tree::anisotropic_refinement_level() const {return _ref_level.copy();}
Array<Int> Tree::coordinates() const {return _coords.copy();}
double Tree::nominal_size() const {return nominal_shape().maxCoeff();}

Mat<> Tree::nominal_shape() const {
  Mat<> nom_shape(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) nom_shape(i_dim) = _root_sz/math::pow(2., _ref_level[i_dim]);
  return nom_shape;
}

Mat<> Tree::nominal_position() const {
  Mat<> pos = _orig;
  auto shape = nominal_shape();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    pos(i_dim) += shape(i_dim)*_coords[i_dim];
  }
  return pos;
}

Mat<> Tree::center() const {return nominal_position() + .5*nominal_shape();}

Tree* Tree::parent() {return _par;}

std::vector<Tree*> Tree::children() {
  std::vector<Tree*> c;
  for (auto& t : _children_storage) c.push_back(t.get());
  return c;
}

std::vector<Tree*> Tree::unique_children() {
  std::vector<Tree*> c;
  for (auto& t : _children_storage) {
    HEXED_ASSERT(t.use_count(), "child is null")
    if (std::none_of(c.begin(), c.end(), [&t](Tree* ptr){return ptr == t.get();})) c.push_back(t.get());
  }
  return c;
}

Tree* Tree::root() {
  Tree* r = this;
  while (!r->is_root()) r = r->parent();
  return r;
}

bool Tree::is_root() {return !_par;}
bool Tree::is_graft() {return root()->_is_graft;}
bool Tree::is_leaf() {return _children_storage.empty();}

bool Tree::is_refined(int i_dim) {
  if (is_leaf()) return false;
  return _children_storage[0] != _children_storage[math::stride(n_dim, 2, i_dim)];
}

bool Tree::has_graft_connection() {
  return std::any_of(_face_connections.begin(), _face_connections.end(), [](_Connection* c)->bool{return c;});
}

std::vector<Tree*> Tree::refine(std::vector<bool> dims) {
  auto ptrs = _refine(dims);
  if (_par) _par->_simplify_aniso_ref();
  return ptrs;
}

std::vector<Tree*> Tree::refine() {
  return refine(std::vector<bool>(n_dim, true));
}

std::vector<Tree*> Tree::refine(int i_dim) {
  std::vector<bool> dims(n_dim, false);
  dims[i_dim] = true;
  return refine(dims);
}

std::vector<Tree*> Tree::unrefine(std::vector<bool> dims) {
  HEXED_ASSERT((int)dims.size() == n_dim, "`refine_dims` has wrong number of entries")
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    HEXED_ASSERT(is_refined(i_dim) || !dims[i_dim], "Cannot unrefine dimension that is not refined.")
    dims[i_dim] = is_refined(i_dim) && !dims[i_dim];
  }
  for (auto& c : _children_storage) HEXED_ASSERT(c->is_leaf(), "At least one child is not a leaf.")
  _children_storage.clear();
  return refine(dims);
}

std::vector<Tree*> Tree::unrefine() {
  return unrefine(std::vector<bool>(n_dim, true));
}

std::vector<Tree*> Tree::unrefine(int i_dim) {
  std::vector<bool> dims(n_dim, false);
  dims[i_dim] = true;
  return unrefine(dims);
}

void Tree::force_unrefine() {_children_storage.clear();}

Tree* Tree::graft(Array<int> ref_level, Array<Int> coords) {
  HEXED_ASSERT(is_root() && !is_graft(), "Can only graft to the root.")
  HEXED_ASSERT(coords.size() == n_dim, "`coords` has wrong number of entries.")
  HEXED_ASSERT(ref_level.size() == n_dim, "`ref_level` has wrong number of entries.")
  _grafts.emplace_back(std::make_unique<Tree>(n_dim, _root_sz, _orig));
  Tree* g = _grafts.back().get();
  g->_ref_level = ref_level.copy();
  g->_coords = coords.copy();
  g->_is_graft = true;
  return g;
}

void Tree::delete_grafts() {
  for (auto& ptr : _grafts) if (ptr) {
    ptr->_clear_connections();
  }
  _clear_connections();
  _grafts.clear();
  _connections.clear();
}

void Tree::connect(std::array<std::vector<Tree*>, 2> trees, Connection_direction dir) {
  HEXED_ASSERT(is_root() && !is_graft(), "Can only add graft connections to the root.")
  for (int i_side = 0; i_side < 2; ++i_side) {
    HEXED_ASSERT((Int)trees[i_side].size() == math::pow(2, n_dim - 1), "Wrong number of trees provided for connection.")
  }
  _connections.emplace_back(std::make_unique<_Connection>(std::array<Tree*, 2>{trees[0][0], trees[1][0]}, dir));
  auto con = _connections.back().get();
  for (int i_side = 0; i_side < 2; ++i_side) {
    auto predicate = [&con, i_side](Tree* t){return t == con->trees[i_side];};
    if (!std::all_of(trees[i_side].begin(), trees[i_side].end(), predicate)) {
      // determine refinement level and coordinates for fake root
      Tree* t0 = trees[i_side][0];
      Array<int> rl = t0->_ref_level.copy();
      Array<Int> coords = t0->_coords.copy();
      for (int i_tree = 0; i_tree < _n_vert()/2; ++i_tree) {
        Tree* t = trees[i_side][i_tree];
        for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim) {
          int j_dim = i_dim + (i_dim >= dir.i_dim[i_side]);
          int coord_diff = math::row_coordinate(n_dim - 1, 2, i_dim, i_tree);
          Tree* neighbor = trees[i_side][i_tree - coord_diff*math::stride(n_dim - 1, 2, i_dim)];
          if (coord_diff && (t != neighbor)) {
            rl[j_dim] = t->_ref_level[j_dim] - 1;
            HEXED_ASSERT(math::mod<Int>(t->_coords[j_dim], 2) == 1,
                         "Coordinate of greater element in refined graft connection must be odd.")
            coords[j_dim] = (t->_coords[j_dim] - 1)/2;
          }
        }
      }
      // check that refinement levels and coordinates of other trees are compatible
      for (int i_tree = 0; i_tree < _n_vert()/2; ++i_tree) {
        Tree* t = trees[i_side][i_tree];
        for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim) {
          int j_dim = i_dim + (i_dim >= dir.i_dim[i_side]);
          int coord_diff = math::row_coordinate(n_dim - 1, 2, i_dim, i_tree);
          Tree* neighbor = trees[i_side][i_tree - math::sign(coord_diff)*math::stride(n_dim - 1, 2, i_dim)];
          int ref = t != neighbor;
          HEXED_ASSERT(t->_ref_level[j_dim] == rl[j_dim] + ref,
                       "Incompatible refinement level in graft connection.")
          HEXED_ASSERT(t->_coords[j_dim] == coords[j_dim] + ref*(coords[j_dim] + coord_diff),
                       "Incompatible coordinates in graft connection.")
        }
      }
      Tree* fake_root = graft(rl, coords);
      con->trees[i_side] = fake_root;
      fake_root->_children_storage.resize(_n_vert());
      for (Row_index ind(n_dim, 2, dir.i_dim[i_side]); ind; ++ind) {
        Tree* t = trees[i_side][ind.i_face_qpoint()];
        HEXED_ASSERT(!t->_par, "Creating a fake parent for a tree that already has one is not yet supported.",
                     assert::Not_implemented_error)
        std::shared_ptr<Tree> child;
        for (Tree* p : t->_fake_parents) if (p) {
          for (std::shared_ptr<Tree>& c : p->_children_storage) if (c.get() == t) child = c;
          HEXED_ASSERT(child.use_count(), "Fake parent/child relationship is not reciprocal.")
        }
        // obtain a shared pointer to `t`
        // without creating any ownership conflicts with existing child or graft pointers
        if (!child.use_count()) {
          for (std::unique_ptr<Tree>& g : _grafts) if (g.get() == t) g.release();
          child.reset(t);
        }
        t->_fake_parents[dir.i_face(i_side)] = fake_root;
        // assign the appropriate children of `fake_root` to point to `t`
        for (int i_row = 0; i_row < 2; ++i_row) fake_root->_children_storage[ind.i_qpoint(i_row)] = child;
      }
    }
    con->trees[i_side]->_face_connections[con->direction.i_face(i_side)] = con;
  }
}

void Tree::connect(std::array<Tree*, 2> trees, Connection_direction dir) {
  std::array<std::vector<Tree*>, 2> vecs;
  int n_tree = math::pow(2, n_dim - 1);
  vecs[0].resize(n_tree, trees[0]);
  vecs[1].resize(n_tree, trees[1]);
  connect(vecs, dir);
}

Tree* Tree::find_leaf(Array<int> ref_level, Array<Int> c, Array<int> b) {
  HEXED_ASSERT(ref_level.size() == n_dim, "`rev_level` has wrong size")
  HEXED_ASSERT(c.size() == n_dim, "`coords` has wrong size")
  HEXED_ASSERT(b.size() >= n_dim, "`bias` has too few entries")
  Array<int> bias = b(0, n_dim).copy();
  // find the relative coordinates in this element's ref level or the specified ref level, whichever is higher
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    int max_level = std::max(_ref_level[i_dim], ref_level[i_dim]);
    Int cell_size = math::pow<Int>(2, max_level - _ref_level[i_dim]);
    Int relative_coord = c[i_dim]*math::pow<Int>(2, max_level - ref_level[i_dim])
                         - bias[i_dim] - _coords[i_dim]*cell_size;
    if (relative_coord < 0 || relative_coord >= cell_size) return nullptr;
  }
  // recursive case: if this element contains the point and has children, one of them should have the element we want
  for (Tree* child : unique_children()) {
    Tree* leaf = child->find_leaf(ref_level, c, bias);
    if (leaf) return leaf;
  }
  // second base case: if the coordinates are in this element, but there are no children,
  // then this is the element we want
  return this;
}

Tree* Tree::find_leaf(int rl, Array<Int> c, Array<int> bias) {
  HEXED_ASSERT(c.size() >= n_dim, "`_coords` has too few elements");
  HEXED_ASSERT(bias.size() >= n_dim, "`bias` has too few elements");
  return find_leaf(Array<int>::make_uniform({n_dim}, rl), c, bias);
}

Tree* Tree::find_leaf(Mat<> nom_pos) {
  HEXED_ASSERT(nom_pos.size() >= n_dim, "`nominal_position` has too few elements");
  Mat<> np = nominal_position();
  Mat<> ns = nominal_shape();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if ((nom_pos(i_dim) < np(i_dim)) || (nom_pos(i_dim) > np(i_dim) + ns(i_dim))) return nullptr;
  }
  for (auto& child : unique_children()) {
    Tree* leaf = child->find_leaf(nom_pos);
    if (leaf) return leaf;
  }
  return this;
}

Tree* Tree::find_neighbor(Array<int> direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  return _neighbor(direction).neighbor;
}

Array<int> Tree::get_direction(int i_face, int n_dim) {
  Array<int> dir = Array<int>::make_uniform({n_dim}, 0);
  dir[i_face/2] = math::sign(i_face%2);
  return dir;
}

Tree* Tree::find_neighbor(int i_face) {
  return find_neighbor(get_direction(i_face, n_dim));
}

// the return value can be 0, 1, or 2, indicating the following:
// 0: there is at least one dimension in which `tree0` has smaller refinement level than `tree1`
// 1: previous condition is false
//    and there is at least one dimension in which `tree1` has smaller refinement level than `tree0`
// 2: the refinement levels are equal
int Tree::_compare_ref_level(Tree* tree0, Tree* tree1, Tree::_Transformation trans) {
  HEXED_ASSERT(tree0 && tree1, "Tree is null.")
  Array<int> rl({2, tree0->n_dim});
  rl(0) = tree0->anisotropic_refinement_level();
  rl(1) = tree1->anisotropic_refinement_level();
  if (trans.used) rl(0) = trans.transform(rl(0));
  Array<int> diff = rl(0) - rl(1);
  if (trans.used) diff[trans.dir.i_dim[!trans.i_side]] = 0; //! \todo make this work for diagonal anisotropic
  if (diff.extreme(0) < 0) return 0;
  if (diff.extreme(1) > 0) return 1;
  return 2;
}

Tree* Tree::_find_parent(int i_face) {
  if (_par) return _par;
  if (i_face >= 0) {
    if (_fake_parents[i_face]) return _fake_parents[i_face];
  }
  HEXED_THROW("No (real or fake) parents.") throw;
}

std::vector<Tree*> Tree::find_neighbors(Array<int> direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  std::vector<Tree*> neighbs;
  // start by finding some leaf neighbor
  auto result = _neighbor(direction);
  if (result.neighbor) {
    int i_face = (direction.abs().sum() == 1) ? result.trans.dir.i_face(!result.trans.i_side) : -1;
    // find a neighbor, not necessarily a leaf, with refinement level not exceeding that of this
    while (_compare_ref_level(this, result.neighbor, result.trans) == 0) {
      result.neighbor = result.neighbor->_find_parent(i_face);
    }
    HEXED_ASSERT(result.neighbor, "Root appears not to satisfy ref level bounds")
    // find all the leaf descendents of that neighbor which are neighbors of this
    Array<int> bias({n_dim});
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      bias[i_dim] = result.direction[i_dim] == 0 ? -1 : result.direction[i_dim] < 0;
    }
    result.neighbor->_add_extremal_levels(neighbs, result.ref_level, result.coords, bias);
    if (neighbs.empty()) neighbs.push_back(result.neighbor); // if neighbor is of lower ref level, add the one neighbor
  }
  return neighbs;
}

std::vector<Tree*> Tree::find_neighbors(int i_face) {
  return find_neighbors(get_direction(i_face, n_dim));
}

Tree::Connection_neighbors Tree::find_connection_neighbors(int i_face) {
  Connection_neighbors neighbors;
  auto result = _neighbor(get_direction(i_face, n_dim));
  HEXED_ASSERT(result.neighbor, "No neighbors on requested face.")
  Tree* search_roots [2] {this, result.neighbor};
  Array<int> rl({2, n_dim});
  Array<Int> coords({2, n_dim});
  for (int i_side = 0; i_side < 2; ++i_side) {
    rl(i_side) = search_roots[i_side]->_ref_level;
    coords(i_side) = search_roots[i_side]->_coords;
  }
  int compare;
  while ((compare = _compare_ref_level(search_roots[0], search_roots[1], result.trans)) != 2) {
    int i_fake_face = result.trans.dir.i_face(result.trans.i_side == compare);
    search_roots[!compare] = search_roots[!compare]->_find_parent(i_fake_face);
  }
  for (bool decrease_rl = true; decrease_rl;) {
    HEXED_ASSERT(rl.extreme(0) >= 0, "negative refinement level")
    // equalize refinement levels so that they are equivalent on both sides of the connection
    for (int i_side = 0; i_side < 2; ++i_side) {
      _Transformation trans = result.trans;
      Array<int> that_rl = rl(!i_side).copy();
      if (trans.used) {
        if (!i_side) trans.reverse();
        that_rl = trans.transform(that_rl);
      }
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) if (i_dim != trans.dir.i_dim[!trans.i_side]) {
        int diff = rl(i_side)[i_dim] - that_rl[i_dim];
        HEXED_ASSERT(diff < 2, "Ref level difference is too large.")
        if (diff == 1) {
          --rl(i_side)[i_dim];
          coords(i_side)[i_dim] -= math::mod<Int>(coords(i_side)[i_dim], 2);
          coords(i_side)[i_dim] /= 2;
        }
      }
    }
    decrease_rl = false;
    // populate neighbors
    for (int i_side = 0; i_side < 2; ++i_side) {
      neighbors.trees[i_side].clear();
      neighbors.trees[i_side].resize(math::pow(2, n_dim - 1), nullptr);
      int j_side = i_side != result.trans.i_side;
      int i_dim = result.trans.dir.i_dim[i_side];
      rl(j_side)[i_dim] = search_roots[j_side]->_ref_level[i_dim];
      coords(j_side)[i_dim] = search_roots[j_side]->_coords[i_dim];
      search_roots[j_side]->_assign_leaves(neighbors.trees[i_side], rl(j_side), coords(j_side),
                                           i_dim, result.trans.dir.face_sign[i_side]);
      // if one of the neighbors is larger than the specified refinement level,
      // we need to try again with a lower refinement level
      for (Tree* n : neighbors.trees[i_side]) if (n) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          if (n->_ref_level[j_dim] < rl(j_side)[j_dim]) {
            decrease_rl = true;
            --rl(j_side)[j_dim];
            coords(j_side)[j_dim] -= math::mod<Int>(coords(j_side)[j_dim], 2);
            coords(j_side)[j_dim] /= 2;
          }
        }
      }
    }
  }
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (Tree* n : neighbors.trees[i_side]) HEXED_ASSERT(n, "Connection neighbor is null.")
  }
  neighbors.direction = result.trans.dir;
  return neighbors;
}

int Tree::count() {
  int total = 1;
  for (auto& child : _children_storage) total += child->count(); //! \todo make this work for aniso
  return total;
}

Array<int> Tree::needs_refine(std::function<bool(Tree*)> include) {
  Array<int> needs({n_dim});
  needs = 0;
  for (int i_face = 0; i_face < 2*n_dim; ++i_face) {
    bool check_con_neighbors = true;
    auto result = _neighbor(get_direction(i_face, n_dim));
    int i_side = result.trans.i_side;
    result.trans.reverse();
    auto neighbors = find_neighbors(i_face);
    auto check_neighbor = [&](Tree* n) {
      if (include(n)) {
        Array<int> desired_rl_diff = result.trans.transform(n->desired_refinement_level()) - desired_refinement_level();
        Array<int> actual_rl_diff = result.trans.transform(n->_ref_level) - _ref_level;
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          needs[i_dim] = needs[i_dim] || (desired_rl_diff[i_dim] > 1 + (i_dim == i_face/2));
          check_con_neighbors = check_con_neighbors && (i_dim == i_face/2 || std::abs(actual_rl_diff[i_dim]) < 2);
        }
      }
    };
    for (Tree* n : neighbors) check_neighbor(n);
    if (check_con_neighbors && !needs.extreme(1) && neighbors.size()) {
      auto con_neighbors = find_connection_neighbors(i_face);
      Array<int> old_needs = needs.copy();
      for (Tree* n : con_neighbors.trees[!i_side]) check_neighbor(n);
      if (!old_needs.equal(needs)) printers::info("(modifying needs ref)", true);
    }
  }
  return needs;
}

void Tree::visualize(std::string format, std::string name) {
  auto vis = Visualizer::create(format, n_dim, n_dim, name, {"tree_level"}, 0., Visualizer::block);
  _visualize(*vis, 0);
}

int Tree::get_status() {
  return _status;
}

void Tree::set_status(int new_status) {
  _status = new_status;
}

void Tree::flood_fill(int new_status) {
  HEXED_ASSERT(new_status != unprocessed, "`flood_fill` may not be used to set _status to `unprocessed`");
  if (!is_leaf()) _children_storage[0]->flood_fill(new_status); // find a leaf element to start
  std::queue<Tree*> to_process;
  to_process.push(this);
  while (!to_process.empty()) {
    Tree* t = to_process.front();
    to_process.pop();
    // if the next element in the queue is already processed, do nothing.
    // if it isn't, set its _status and add all it's neighbors to the queue.
    // this could result in some cells being in the queue multiple times,
    // but that's not a problem.
    if (t->_status == unprocessed) {
      t->_status = new_status;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int sign : {-1, 1}) {
          Array<int> direct({n_dim});
          direct = 0;
          direct[i_dim] = sign;
          for (Tree* neighb : t->find_neighbors(direct)) {
            to_process.push(neighb);
          }
        }
      }
    }
  }
}

void Tree::clear_status() {
  _status = unprocessed;
  for (auto& child : _children_storage) child->clear_status();
}

void Tree::_add_extremal_levels(std::vector<Tree*>& add_to, Array<int> ref_level, Array<Int> coords, Array<int> bias) {
  if (is_leaf()) add_to.push_back(this);
  else {
    std::vector<Tree*> added;
    for (int i_child = 0; i_child < math::pow(2, n_dim); ++i_child) {
      Tree* child = _children_storage[i_child].get();
      bool add = true;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        Int scale0 = math::pow<Int>(2, ref_level[i_dim]);
        Int scale1 = math::pow<Int>(2, child->_ref_level[i_dim]);
        bool overlap = scale0*child->_coords[i_dim] < scale1*(coords[i_dim] + 1)
                       && scale0*(child->_coords[i_dim] + 1) > scale1*coords[i_dim];
        add = add && (bias[i_dim] == -1 || bias[i_dim] == math::row_coordinate(n_dim, 2, i_dim, i_child))
                  && (bias[i_dim] != -1 || overlap);
      }
      add = add && std::none_of(added.begin(), added.end(), [child](Tree* t){return t == child;});
      if (add) {
        added.push_back(child);
        child->_add_extremal_levels(add_to, ref_level, coords, bias);
      }
    }
  }
}

void Tree::_assign_leaves(std::vector<Tree*>& assign_to, Array<int> ref_level, Array<Int> coords, int i_dim, int sign) {
  if (!is_leaf()) {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      bool any_ref = false;
      bool all_ref = true;
      for (Row_index index(n_dim, 2, j_dim); index; ++index) {
        bool ref = _children_storage[index.i_qpoint(0)].get() != _children_storage[index.i_qpoint(1)].get();
        any_ref = any_ref || ref;
        all_ref = all_ref && ref;
      }
      HEXED_ASSERT(any_ref == all_ref, "refinement mismatch")
    }
  }
  for (Row_index index(n_dim, 2, i_dim); index; ++index) {
    int i_child = index.i_qpoint(sign);
    if (is_leaf()) {
      bool assign = true;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        int row = math::row_coordinate(n_dim, 2, j_dim, i_child);
        Int scale0 = math::pow<Int>(2, ref_level[j_dim]);
        Int scale1 = math::pow<Int>(2, _ref_level[j_dim]);
        assign = assign && ((_coords[j_dim] + row)*scale0 == (coords[j_dim] + row)*scale1);
      }
      if (assign) assign_to[index.i_face_qpoint()] = this;
    } else {
      _children_storage[i_child]->_assign_leaves(assign_to, ref_level, coords, i_dim, sign);
    }
  }
}

std::vector<Tree*> Tree::_refine(std::vector<bool> dims) {
  HEXED_ASSERT(is_leaf(), "can only refine leaf")
  HEXED_ASSERT((Int)dims.size() == n_dim, "`refine_dims` has wrong number of entries")
  int n_child = math::pow(2, n_dim);
  if (std::none_of(dims.begin(), dims.end(), [](bool b){return b;})) return std::vector<Tree*>(n_child, this);
  _children_storage.resize(n_child);
  for (int i_child = 0; i_child < n_child; ++i_child) {
    bool redundant = false;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      redundant = redundant || (!dims[i_dim] && math::row_coordinate(n_dim, 2, i_dim, i_child));
    }
    if (redundant) continue;
    auto child = std::make_shared<Tree>(n_dim, _root_sz, _orig);
    child->_par = this;
    child->_ref_level = _ref_level;
    child->_coords = _coords;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) if (dims[i_dim]) {
      child->_ref_level[i_dim] += 1;
      child->_coords[i_dim] *= 2;
      child->_coords[i_dim] += math::row_coordinate(n_dim, 2, i_dim, i_child);
    }
    for (int j_child = 0; j_child < n_child; ++j_child) {
      bool assign = true;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        int diff = math::row_coordinate(n_dim, 2, i_dim, j_child) - math::row_coordinate(n_dim, 2, i_dim, i_child);
        assign = assign && !(dims[i_dim] && diff);
      }
      if (assign) _children_storage[j_child] = child;
    }
  }
  return children();
}

void Tree::_interchange_aniso_ref() {
  if (is_leaf()) return;
  int n_child = math::pow(2, n_dim);
  std::vector<bool> any_refined(n_dim, false);
  std::vector<bool> all_refined(n_dim, true);
  bool allowed = true;
  for (auto& child : _children_storage) {
    for (_Connection* c : child->_face_connections) allowed = allowed && !c;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      bool ref = child->is_refined(i_dim);
      any_refined[i_dim] = any_refined[i_dim] || ref;
      all_refined[i_dim] = all_refined[i_dim] && ref;
    }
  }
  if (!allowed) return;
  std::vector<bool> pass_down(n_dim);
  std::vector<bool> retain(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if (any_refined[i_dim] != all_refined[i_dim]) return;
    pass_down[i_dim] = is_refined(i_dim) && !any_refined[i_dim];
    retain[i_dim] = is_refined(i_dim) && all_refined[i_dim];
  }
  auto b = [](bool arg){return arg;};
  if (std::none_of(pass_down.begin(), pass_down.end(), b) || std::none_of(retain.begin(), retain.end(), b)) return;
  std::vector<std::vector<std::shared_ptr<Tree>>> grandchildren;
  for (auto& child : _children_storage) {
    grandchildren.push_back(child->_children_storage);
  }
  force_unrefine();
  _refine(retain);
  for (int i_child = 0; i_child < n_child; ++i_child) {
    auto child = _children_storage[i_child];
    child->_children_storage.resize(n_child); // does nothing if `child` has already been visited
    for (int j_child = 0; j_child < n_child; ++j_child) {
      bool assign = true;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        int coord = math::row_coordinate(n_dim, 2, i_dim, i_child);
        assign = assign && (is_refined(i_dim) || coord == math::row_coordinate(n_dim, 2, i_dim, j_child));
      }
      if (assign) {
        child->_children_storage[j_child] = grandchildren[i_child][j_child];
        grandchildren[i_child][j_child]->_par = child.get();
      }
    }
  }
}

void Tree::_collapse_aniso_ref() {
  if (is_leaf()) return;
  HEXED_ASSERT(_children_storage.size() > 1, "one child")
  bool collapse = true;
  while (collapse) {
    for (auto& child : _children_storage) {
      collapse = collapse && !child->is_leaf();
      for (_Connection* c : child->_face_connections) collapse = collapse && !c;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        collapse = collapse && !(child->is_refined(i_dim) && is_refined(i_dim));
      }
    }
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      bool any = false;
      bool all = true;
      for (auto& child : _children_storage) {
        any = any || child->is_refined(i_dim);
        all = all && child->is_refined(i_dim);
      }
      collapse = collapse && (any == all);
    }
    if (collapse) {
      for (int i_child = 0; i_child < (int)_children_storage.size(); ++i_child) {
        _children_storage[i_child] = _children_storage[i_child]->_children_storage[i_child];
        _children_storage[i_child]->_par = this;
      }
    }
  }
}

void Tree::_simplify_aniso_ref() {
  _collapse_aniso_ref();
  _interchange_aniso_ref();
}

Array<int> Tree::_Transformation::transform(Array<int> ref_level, bool rot) {
  int i_dim = dir.i_dim[!i_side];
  int j_dim = dir.i_dim[i_side];
  Array<int> transformed = ref_level - this_root->_ref_level;
  std::swap(transformed[i_dim], transformed[j_dim]);
  int dim0 = i_dim == 0;
  int dim1 = 1 + (i_dim <= 1);
  if (rot && dir.rotate%2) std::swap(transformed[dim0], transformed[dim1]);
  transformed += that_root->_ref_level;
  return transformed;
}

void Tree::_Transformation::reverse() {
  i_side = !i_side;
  std::swap(this_root, that_root);
}

Tree::_Neighbor_result Tree::_neighbor(Array<int> dir_arg) {
  Array<int> direction = dir_arg.copy();
  // compute the coordinates and bias which will identify the neighbor
  Array<int> bias({n_dim});
  Array<Int> coords({n_dim});
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    coords[i_dim] = _coords[i_dim] + (direction[i_dim] > 0);
    bias[i_dim] = (direction[i_dim] < 0);
  }
  Tree* r = root();
  Tree* n = nullptr;
  int i_face = -1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if (direction[i_dim]) {
      if (i_face == -1) i_face = 2*i_dim + (direction[i_dim] > 0);
      else i_face = -2;
    }
  }
  _Transformation trans {
    .used = i_face >= 0,
    .this_root = r,
    .that_root = r,
    .dir{{i_face/2, i_face/2}, {1, 0}},
    .i_side = !(i_face%2),
  };
  Array<int> ref_level = _ref_level.copy();
  Array<Int> search_coords = _coords.copy();
  if (i_face >= 0) if (!n) {
    Tree* search_root = this;
    Array<int> search_bias({n_dim});
    Array<int> search_direction({n_dim});
    while (search_root) {
      if (search_root->_fake_parents[i_face]) {
        search_root = search_root->_fake_parents[i_face];
      } else if (search_root->_face_connections[i_face]) {
        Int scale = math::pow<Int>(2, _ref_level[i_face/2] - search_root->_ref_level[i_face/2]);
        if (_coords[i_face/2] + i_face%2 != (search_root->_coords[i_face/2] + i_face%2)*scale) {
          search_root = nullptr;
          break;
        }
        Tree* this_root = search_root;
        auto trees = this_root->_face_connections[i_face]->trees;
        trans.i_side = trees[1] == this_root;
        search_root = trees[!trans.i_side];
        trans.this_root = this_root;
        trans.that_root = search_root;
        search_coords = coords;
        for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
          int rl_diff = _ref_level[k_dim] - this_root->_ref_level[k_dim];
          search_coords[k_dim] -= this_root->_coords[k_dim]*math::pow<Int>(2, rl_diff);
        }
        trans.dir = this_root->_face_connections[i_face]->direction;
        int i_dim = trans.dir.i_dim[!trans.i_side];
        int j_dim = trans.dir.i_dim[trans.i_side];
        search_coords[j_dim] = search_coords[i_dim];
        search_coords[i_dim] = trans.dir.face_sign[!trans.i_side];
        ref_level = _ref_level - this_root->_ref_level;
        ref_level[j_dim] = ref_level[i_dim];
        ref_level[i_dim] = 0;
        search_bias = 0;
        search_bias[i_dim] = trans.dir.face_sign[!trans.i_side];
        search_direction = direction;
        search_direction[j_dim] = 0;
        search_direction[i_dim] = math::sign(!trans.dir.face_sign[!trans.i_side]);
        if (trans.dir.flip_tangential()) { // implies different dims
          search_coords[j_dim] = math::pow<Int>(2, ref_level[j_dim]) - search_coords[j_dim] - 1;
        }
        if (n_dim == 3) {
          int dim0 = (i_dim + 1)%3;
          int dim1 = (i_dim + 2)%3;
          int rotation_sign = math::sign(trans.dir.face_sign[!trans.i_side] != trans.dir.face_sign[0]);
          int n_rot = math::mod(trans.dir.rotate*rotation_sign, 4);
          for (int rot = 0; rot < n_rot; ++rot) {
            // achieve rotation by transposing dimensions and then inverting one of them
            // first transpose
            std::swap(search_coords[dim0], search_coords[dim1]);
            std::swap(search_bias[dim0], search_bias[dim1]);
            std::swap(ref_level[dim0], ref_level[dim1]);
            // now flip
            search_coords[dim1] = math::pow<Int>(2, ref_level[dim1]) - search_coords[dim1] - 1;
          }
        }
        for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
          search_coords[k_dim] += search_root->_coords[k_dim]*math::pow<Int>(2, ref_level[k_dim]);
        }
        ref_level += search_root->_ref_level;
        break;
      } else {
        search_root = search_root->_par;
      }
    }
    if (search_root) {
      n = search_root->find_leaf(ref_level, search_coords, search_bias);
      if (n) {
        direction = search_direction;
        trans.used = true;
      }
    }
  }
  if (!n) n = r->find_leaf(_ref_level, coords, bias);
  return {n, direction, trans, ref_level, search_coords};
}

void Tree::_clear_connections() {
  for (_Connection*& c : _face_connections) c = nullptr;
  for (Tree* t : unique_children()) t->_clear_connections();
}

void Tree::_visualize(Visualizer& vis, int tree_level) {
  int nv = math::pow(2, n_dim);
  std::vector<Int> shape(n_dim, 2);
  shape.insert(shape.begin(), 1);
  Array<double> data = Array<double>::make_uniform(shape, tree_level);
  shape[0] = n_dim;
  Array<double> pos(shape);
  for (int i_point = 0; i_point < nv; ++i_point) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      int rc = math::row_coordinate(n_dim, 2, i_dim, i_point);
      pos(i_dim)[i_point] = nominal_position()[i_dim] + rc*nominal_shape()[i_dim];
    }
  }
  vis.write_block(pos, data);
  for (Tree* child : unique_children()) child->_visualize(vis, tree_level + 1);
}

}
