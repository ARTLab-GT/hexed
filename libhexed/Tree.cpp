#include <queue>
#include <hexed/Tree.hpp>
#include <hexed/Row_index.hpp>

namespace hexed {

std::array<std::vector<Element*>, 2> Tree::Connection_neighbors::elements() {
  std::array<std::vector<Element*>, 2> elems;
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (Tree* t : trees[i_side]) elems[i_side].push_back(t->elem.get());
  }
  return elems;
}

Tree::Tree(int nd, double root_size, Mat<> origin)
: n_dim{nd}
, elem(this)
, _root_sz{root_size}
, _ref_level{Array<int>::make_uniform({nd}, 0)}, _coords{Eigen::VectorXi::Zero(nd)}
, _par{nullptr}
, _children_storage()
, _face_connections(2*n_dim, nullptr)
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
Array<int> Tree::anisotropic_refinement_level() const {return _ref_level.copy();}
Eigen::VectorXi Tree::coordinates() const {return _coords;}
double Tree::nominal_size() const {return nominal_shape().maxCoeff();}

Mat<> Tree::nominal_shape() const {
  Mat<> nom_shape(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) nom_shape(i_dim) = _root_sz/math::pow(2, _ref_level[i_dim]);
  return nom_shape;
}

Mat<> Tree::nominal_position() const {return nominal_shape().cwiseProduct(_coords.cast<double>()) + _orig;}
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

void Tree::unrefine(std::vector<bool> dims) {
  HEXED_ASSERT((int)dims.size() == n_dim, "`refine_dims` has wrong number of entries")
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    HEXED_ASSERT(is_refined(i_dim) || !dims[i_dim], "Cannot unrefine dimension that is not refined.")
    dims[i_dim] = is_refined(i_dim) && !dims[i_dim];
  }
  for (auto& c : _children_storage) HEXED_ASSERT(c->is_leaf(), "At least one child is not a leaf.")
  _children_storage.clear();
  refine(dims);
}

void Tree::unrefine() {
  unrefine(std::vector<bool>(n_dim, true));
}

void Tree::unrefine(int i_dim) {
  std::vector<bool> dims(n_dim, false);
  dims[i_dim] = true;
  unrefine(dims);
}

void Tree::force_unrefine() {_children_storage.clear();}

Tree* Tree::graft(Array<int> ref_level, Eigen::VectorXi coords) {
  HEXED_ASSERT(is_root(), "Can only graft to the root.")
  HEXED_ASSERT(coords.size() == n_dim, "`coords` has wrong number of entries.")
  HEXED_ASSERT(ref_level.size() == n_dim, "`ref_level` has wrong number of entries.")
  _grafts.emplace_back(std::make_unique<Tree>(n_dim, _root_sz, _orig));
  Tree* g = _grafts.back().get();
  g->_ref_level = ref_level.copy();
  g->_coords = coords;
  g->_is_graft = true;
  return g;
}

void Tree::delete_grafts() {
  for (auto& ptr : _grafts) ptr->_clear_connections();
  _grafts.clear();
  _connections.clear();
  _clear_connections();
}

void Tree::connect(std::array<std::vector<Tree*>, 2> trees, Connection_direction dir) {
  HEXED_ASSERT(is_root(), "Can only add graft connections to the root.")
  for (int i_side = 0; i_side < 2; ++i_side) {
    HEXED_ASSERT((Int)trees[i_side].size() == math::pow(2, n_dim - 1), "Wrong number of trees provided for connection.")
  }
  _connections.emplace_back(std::make_unique<_Connection>(std::array<Tree*, 2>{trees[0][0], trees[1][0]}, dir));
  auto con = _connections.back().get();
  for (int i_side = 0; i_side < 2; ++i_side) {
    auto predicate = [&con, i_side](Tree* t){return t == con->trees[i_side];};
    if (std::all_of(trees[i_side].begin(), trees[i_side].end(), predicate)) {
    } else {
      HEXED_THROW("Trees must all be the same.")
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

void delete_grafts() {
}

Tree* Tree::find_leaf(Array<int> ref_level, Eigen::VectorXi c, Eigen::VectorXi b) {
  HEXED_ASSERT(ref_level.size() == n_dim, "`rev_level` has wrong size")
  HEXED_ASSERT(c.size() == n_dim, "`coords` has wrong size")
  HEXED_ASSERT(b.size() >= n_dim, "`bias` has too few entries")
  Eigen::VectorXi bias = b(Eigen::seqN(0, n_dim));
  // find the relative coordinates in this element's ref level or the specified ref level, whichever is higher
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    int max_level = std::max(_ref_level[i_dim], ref_level[i_dim]);
    int cell_size = math::pow(2, max_level - _ref_level[i_dim]);
    int relative_coord = c(i_dim)*math::pow(2, max_level - ref_level[i_dim]) - bias(i_dim) - _coords(i_dim)*cell_size;
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

Tree* Tree::find_leaf(int rl, Eigen::VectorXi c, Eigen::VectorXi bias) {
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

Tree* Tree::find_neighbor(Eigen::VectorXi direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  return _neighbor(direction).neighbor;
}

Eigen::VectorXi Tree::get_direction(int i_face, int n_dim) {
  Eigen::VectorXi dir = Eigen::VectorXi::Zero(n_dim);
  dir(i_face/2) = math::sign(i_face%2);
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

std::vector<Tree*> Tree::find_neighbors(Eigen::VectorXi direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  std::vector<Tree*> neighbs;
  // start by finding some leaf neighbor
  auto result = _neighbor(direction);
  if (result.neighbor) {
    // find a neighbor, not necessarily a leaf, with refinement level not exceeding that of this
    while (_compare_ref_level(this, result.neighbor, result.trans) == 0) result.neighbor = result.neighbor->parent();
    HEXED_ASSERT(result.neighbor, "Root appears not to satisfy ref level bounds")
    // find all the leaf descendents of that neighbor which are neighbors of this
    Eigen::VectorXi bias(n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      bias(i_dim) = result.direction(i_dim) == 0 ? -1 : result.direction(i_dim) < 0;
    }
    result.neighbor->_add_extremal_levels(neighbs, bias);
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
  Tree* search_roots [2];
  search_roots[0] = this;
  search_roots[1] = result.neighbor;
  int compare;
  while ((compare = _compare_ref_level(search_roots[0], search_roots[1], result.trans)) != 2) {
    search_roots[!compare] = search_roots[!compare]->_par;
  }
  for (int i_side = 0; i_side < 2; ++i_side) {
    neighbors.trees[i_side].resize(math::pow(2, n_dim - 1), nullptr);
    int j_side = i_side != result.trans.i_side;
    search_roots[j_side]->_assign_leaves(neighbors.trees[i_side], search_roots[j_side],
                                         result.trans.dir.i_dim[i_side], result.trans.dir.face_sign[i_side]);
    for (Tree* n : neighbors.trees[i_side]) HEXED_ASSERT(n, "null element returned")
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
    auto result = _neighbor(get_direction(i_face, n_dim));
    result.trans.reverse();
    auto neighbors = find_neighbors(i_face);
    for (Tree* n : neighbors) if (include(n)) {
      Array<int> rl_diff = result.trans.transform(n->_ref_level) - _ref_level;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        needs[i_dim] = needs[i_dim] || (rl_diff[i_dim] > 1);
      }
    }
  }
  return needs;
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
          Eigen::VectorXi direct(n_dim);
          direct.setZero();
          direct(i_dim) = sign;
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

void Tree::_add_extremal_levels(std::vector<Tree*>& add_to, Eigen::VectorXi bias) {
  if (is_leaf()) add_to.push_back(this);
  else {
    std::vector<Tree*> added;
    for (int i_child = 0; i_child < math::pow(2, n_dim); ++i_child) {
      bool add = true;
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        add = add && (bias(i_dim) == -1 || bias(i_dim) == math::row_coordinate(n_dim, 2, i_dim, i_child));
      }
      Tree* child = _children_storage[i_child].get();
      add = add && std::none_of(added.begin(), added.end(), [child](Tree* t){return t == child;});
      if (add) {
        added.push_back(child);
        child->_add_extremal_levels(add_to, bias);
      }
    }
  }
}

void Tree::_assign_leaves(std::vector<Tree*>& assign_to, Tree* search_root, int i_dim, int sign) {
  for (Row_index index(n_dim, 2, i_dim); index; ++index) {
    int i_child = index.i_qpoint(sign);
    if (is_leaf()) {
      bool assign = true;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        int row = math::row_coordinate(n_dim, 2, j_dim, i_child);
        int scale = math::pow(2, _ref_level[j_dim] - search_root->_ref_level[j_dim]);
        assign = assign && (_coords(j_dim) + row == (search_root->_coords(j_dim) + row)*scale);
      }
      if (assign) assign_to[index.i_face_qpoint()] = this;
    } else {
      _children_storage[i_child]->_assign_leaves(assign_to, search_root, i_dim, sign);
    }
  }
}

std::vector<Tree*> Tree::_refine(std::vector<bool> dims) {
  HEXED_ASSERT(is_leaf(), "can only refine leaf")
  HEXED_ASSERT((int)dims.size() == n_dim, "`refine_dims` has wrong number of entries")
  if (std::none_of(dims.begin(), dims.end(), [](bool b){return b;})) return {};
  int n_child = math::pow(2, n_dim);
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
  return unique_children();
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

Array<int> Tree::_Transformation::transform(Array<int> ref_level) {
  int i_dim = dir.i_dim[!i_side];
  int j_dim = dir.i_dim[i_side];
  Array<int> transformed = ref_level + that_root->_ref_level - this_root->_ref_level;
  transformed[j_dim] = ref_level[i_dim] + that_root->_ref_level[j_dim]
                                        - this_root->_ref_level[i_dim];
  transformed[i_dim] = ref_level[j_dim] + that_root->_ref_level[i_dim]
                                        - this_root->_ref_level[j_dim];
  return transformed;
}

void Tree::_Transformation::reverse() {
  i_side = !i_side;
  std::swap(this_root, that_root);
}

Tree::_Neighbor_result Tree::_neighbor(Eigen::VectorXi direction) {
  // compute the coordinates and bias which will identify the neighbor
  Eigen::VectorXi bias(n_dim);
  Eigen::VectorXi coords(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    coords(i_dim) = _coords[i_dim] + (direction(i_dim) > 0);
    bias(i_dim) = (direction(i_dim) < 0);
  }
  Tree* r = root();
  Tree* n = nullptr;
  int i_face = -1;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if (direction(i_dim)) {
      if (i_face == -1) i_face = 2*i_dim + (direction(i_dim) > 0);
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
  if (i_face >= 0) if (!n) {
    Tree* search_root = this;
    Array<int> ref_level({n_dim});
    Eigen::VectorXi search_coords(n_dim);
    Eigen::VectorXi search_bias(n_dim);
    Eigen::VectorXi search_direction(n_dim);
    while (search_root) {
      if (search_root->_face_connections[i_face]) {
        int scale = math::pow(2, _ref_level[i_face/2] - search_root->_ref_level[i_face/2]);
        if (_coords(i_face/2) + i_face%2 != (search_root->_coords(i_face/2) + i_face%2)*scale) {
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
          search_coords(k_dim) -= this_root->_coords(k_dim)*math::pow(2, rl_diff);
        }
        trans.dir = this_root->_face_connections[i_face]->direction;
        int i_dim = trans.dir.i_dim[!trans.i_side];
        int j_dim = trans.dir.i_dim[trans.i_side];
        search_coords(j_dim) = search_coords(i_dim);
        search_coords(i_dim) = trans.dir.face_sign[!trans.i_side];
        ref_level = trans.transform(_ref_level);
        ref_level[i_dim] = search_root->_ref_level[i_dim];
        search_bias.setZero();
        search_bias(i_dim) = trans.dir.face_sign[!trans.i_side];
        search_direction = direction;
        search_direction(j_dim) = 0;
        search_direction(i_dim) = math::sign(!trans.dir.face_sign[!trans.i_side]);
        if (trans.dir.flip_tangential()) { // implies different dims
          int rl_diff = ref_level[j_dim] - search_root->_ref_level[j_dim];
          search_coords(j_dim) = math::pow(2, rl_diff) - search_coords(j_dim);
          search_bias(j_dim) = 1;
        }
        for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
          int rl_diff = ref_level[k_dim] - search_root->_ref_level[k_dim];
          search_coords(k_dim) += search_root->_coords(k_dim)*math::pow(2, rl_diff);
        }
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
  return {n, direction, trans};
}

void Tree::_clear_connections() {
  for (_Connection*& c : _face_connections) c = nullptr;
  for (Tree* t : unique_children()) t->_clear_connections();
}

}
