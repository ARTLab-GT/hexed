#include <Tree.hpp>
#include <queue>

namespace hexed {

void Tree::_add_extremal_levels(std::vector<Tree*>& add_to, Eigen::VectorXi bias) {
  if (is_leaf()) add_to.push_back(this);
  else for (auto& child : _children_storage) {
    bool add = true;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      add = add && (bias(i_dim) == -1 || bias(i_dim) == child->_coords(i_dim)%2);
    }
    if (add) child->_add_extremal_levels(add_to, bias);
  }
}

void Tree::_interchange_aniso_ref() {
  if (is_leaf()) return;
  std::vector<bool> any_refined(n_dim, false);
  std::vector<bool> all_refined(n_dim, true);
  for (auto& child : _children_storage) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      bool ref = child->is_refined(i_dim);
      any_refined[i_dim] = any_refined[i_dim] || ref;
      all_refined[i_dim] = all_refined[i_dim] && ref;
    }
  }
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
  unrefine();
  refine(retain);
  int n_child = grandchildren.size();
  for (int i_child = 0; i_child < n_child; ++i_child) {
    auto child = _children_storage[i_child];
    child->_children_storage.resize(n_child);
    for (int j_child = 0; j_child < n_child; ++j_child) {
      if (i_child == j_child || _children_storage[i_child] != _children_storage[j_child]) {
        child->_children_storage[i_child] = grandchildren[i_child][j_child];
        grandchildren[i_child][j_child]->_par = child.get();
      }
    }
  }
}

void Tree::_collapse_aniso_ref() {
  if (is_leaf()) return;
  bool collapse = true;
  while (collapse) {
    for (auto& child : _children_storage) {
      collapse = collapse && !child->is_leaf();
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        collapse = collapse && (child->is_refined(i_dim) != is_refined(i_dim));
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
  _interchange_aniso_ref();
  _collapse_aniso_ref();
}

Tree::Tree(int nd, double root_size, Mat<> origin)
: _root_sz{root_size}, _ref_level{Array<int>::make_uniform({nd}, 0)}, _coords{Eigen::VectorXi::Zero(nd)}
, _par{nullptr}, _children_storage()
, _status{unprocessed}
, n_dim{nd}
, elem(this)
{
  HEXED_ASSERT(origin.size() >= n_dim, "`origin` is too small");
  _orig = origin(Eigen::seqN(0, n_dim));
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

bool Tree::is_root() const {return !_par;}
bool Tree::is_leaf() const {return _children_storage.empty();}

bool Tree::is_refined(int i_dim) const {
  if (is_leaf()) return false;
  return _children_storage[0] != _children_storage[math::stride(n_dim, 2, i_dim)];
}

void Tree::refine(std::vector<bool> dims) {
  HEXED_ASSERT(is_leaf(), "can only refine leaf")
  std::cout << "refining: ";
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
      if (assign) std::cout << format_str("(%i %i)", i_child, j_child);
    }
  }
  std::cout << std::endl;
  if (_par) _par->_simplify_aniso_ref();
}

void Tree::refine() {
  refine(std::vector<bool>(n_dim, true));
}

void Tree::refine(int i_dim) {
  std::vector<bool> dims(n_dim, false);
  dims[i_dim] = true;
  refine(dims);
}

void Tree::unrefine() {_children_storage.clear();}

Tree* Tree::find_leaf(int rl, Eigen::VectorXi c, Eigen::VectorXi bias) {
  HEXED_ASSERT(c.size() >= n_dim, "`_coords` has too few elements");
  HEXED_ASSERT(bias.size() >= n_dim, "`bias` has too few elements");
  // find the relative coordinates in this element's ref level or the specified ref level, whichever is higher
  c = c(Eigen::seqN(0, n_dim));
  int max_level = std::max(refinement_level(), rl);
  int cell_size = math::pow(2, max_level - refinement_level());
  Eigen::MatrixXi relative_coords = c*math::pow(2, max_level - rl) - bias(Eigen::seqN(0, n_dim)) - _coords*cell_size;
  // first base case: if the coordinates are outside this element, return null
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if ((relative_coords(i_dim) < 0) || (relative_coords(i_dim) >= cell_size)) return nullptr;
  }
  // recursive case: if this element contains the point and has children, one of them should have the element we want
  for (auto& child : _children_storage) {
    Tree* leaf = child->find_leaf(rl, c, bias);
    if (leaf) return leaf;
  }
  // second base case: if the coordinates are in this element, but there are no children,
  // then this is the element we want
  return this;
}

Tree* Tree::find_leaf(Mat<> nom_pos) {
  HEXED_ASSERT(nom_pos.size() >= n_dim, "`nominal_position` has too few elements");
  Mat<> np = nominal_position();
  double ns = nominal_size();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if ((nom_pos(i_dim) < np(i_dim)) || (nom_pos(i_dim) > np(i_dim) + ns)) return nullptr;
  }
  for (auto& child : _children_storage) {
    Tree* leaf = child->find_leaf(nom_pos);
    if (leaf) return leaf;
  }
  return this;
}

Tree* Tree::find_neighbor(Eigen::VectorXi direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  // compute the coordinates and bias which will identify the neighbor
  Eigen::VectorXi bias(n_dim);
  Eigen::VectorXi c(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    c(i_dim) = _coords[i_dim] + (direction(i_dim) > 0);
    bias(i_dim) = (direction(i_dim) < 0);
  }
  // use `find_leaf` on the root element to find the neighbor
  return root()->find_leaf(refinement_level(), c, bias);
}

std::vector<Tree*> Tree::find_neighbors(Eigen::VectorXi direction) {
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  std::vector<Tree*> neighbs;
  // start by finding some leaf neighbor
  Tree* main_neighbor = find_neighbor(direction);
  if (main_neighbor) {
    // find a neighbor, not necessarily a leaf, with the same refinement level as this
    while (main_neighbor->refinement_level() > refinement_level()) main_neighbor = main_neighbor->parent();
    // find all the leaf descendents of that neighbor which are neighbors of this
    Eigen::VectorXi bias(n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) bias(i_dim) = direction(i_dim) == 0 ? -1 : direction(i_dim) < 0;
    main_neighbor->_add_extremal_levels(neighbs, bias);
  }
  return neighbs;
}

int Tree::count() {
  int total = 1;
  for (auto& child : _children_storage) total += child->count();
  return total;
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

}
