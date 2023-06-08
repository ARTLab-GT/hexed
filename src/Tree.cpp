#include <Tree.hpp>

namespace hexed
{

void Tree::add_extremal_leves(std::vector<Tree*>& add_to, Eigen::VectorXi bias)
{
  if (is_leaf()) add_to.push_back(this);
  else for (auto& child : children_storage) {
    bool add = true;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      add = add && (bias(i_dim) == -1 || bias(i_dim) == child->coords(i_dim)%2);
    }
    if (add) child->add_extremal_leves(add_to, bias);
  }
}

Tree::Tree(int nd, double root_size, Mat<> origin)
: root_sz{root_size}, ref_level{0}, coords{Eigen::VectorXi::Zero(nd)},
  par{nullptr}, children_storage(),
  n_dim{nd}
{
  HEXED_ASSERT(origin.size() >= n_dim, "`origin` is too small");
  orig = origin(Eigen::seqN(0, n_dim));
}

Mat<> Tree::origin() const {return orig;}
int Tree::refinement_level() const {return ref_level;}
Eigen::VectorXi Tree::coordinates() const {return coords;}
double Tree::nominal_size() const {return root_sz/math::pow(2, ref_level);}
Mat<> Tree::nominal_position() const {return nominal_size()*coords.cast<double>() + orig;}

Tree* Tree::parent() {return par;}

std::vector<Tree*> Tree::children()
{
  std::vector<Tree*> c;
  for (auto& t : children_storage) c.push_back(t.get());
  return c;
}

Tree* Tree::root()
{
  Tree* r = this;
  while (!r->is_root()) r = r->parent();
  return r;
}

bool Tree::is_root() const {return !par;}
bool Tree::is_leaf() const {return children_storage.empty();}

void Tree::refine()
{
  for (int i_child = 0; i_child < math::pow(2, n_dim); ++i_child) {
    children_storage.emplace_back(new Tree(n_dim, root_sz, orig));
    Tree& child = *children_storage.back();
    child.par = this;
    child.ref_level = ref_level + 1;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      int stride = math::pow(2, n_dim - 1 - i_dim);
      child.coords(i_dim) = 2*coords[i_dim] + (i_child/stride)%2;
    }
  }
}

void Tree::unrefine() {children_storage.clear();}

Tree* Tree::find_leaf(int rl, Eigen::VectorXi c, Eigen::VectorXi bias)
{
  HEXED_ASSERT(c.size() >= n_dim, "`coords` has too few elements");
  HEXED_ASSERT(bias.size() >= n_dim, "`bias` has too few elements");
  // find the relative coordinates in this element's ref level or the specified ref level, whichever is higher
  c = c(Eigen::seqN(0, n_dim));
  int max_level = std::max(ref_level, rl);
  int cell_size = math::pow(2, max_level - ref_level);
  Eigen::MatrixXi relative_coords = c*math::pow(2, max_level - rl) - bias(Eigen::seqN(0, n_dim)) - coords*cell_size;
  // first base case: if the coordinates are outside this element, return null
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if ((relative_coords(i_dim) < 0) || (relative_coords(i_dim) >= cell_size)) return nullptr;
  }
  // recursive case: if this element contains the point and has children, one of them should have the element we want
  for (auto& child : children_storage) {
    Tree* leaf = child->find_leaf(rl, c, bias);
    if (leaf) return leaf;
  }
  // second base case: if the coordinates are in this element, but there are no children, then this is the element we want
  return this;
}

Tree* Tree::find_leaf(Mat<> nom_pos)
{
  HEXED_ASSERT(nom_pos.size() >= n_dim, "`nominal_position` has too few elements");
  Mat<> np = nominal_position();
  double ns = nominal_size();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    if ((nom_pos(i_dim) < np(i_dim)) || (nom_pos(i_dim) > np(i_dim) + ns)) return nullptr;
  }
  for (auto& child : children_storage) {
    Tree* leaf = child->find_leaf(nom_pos);
    if (leaf) return leaf;
  }
  return this;
}

Tree* Tree::find_neighbor(Eigen::VectorXi direction)
{
  HEXED_ASSERT(direction.size() >= n_dim, "`direction` has too few elements");
  // compute the coordinates and bias which will identify the neighbor
  Eigen::VectorXi bias(n_dim);
  Eigen::VectorXi c(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    c(i_dim) = coords[i_dim] + (direction(i_dim) > 0);
    bias(i_dim) = (direction(i_dim) < 0);
  }
  // use `find_leaf` on the root element to find the neighbor
  return root()->find_leaf(ref_level, c, bias);
}

std::vector<Tree*> Tree::find_neighbors(Eigen::VectorXi direction)
{
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
    main_neighbor->add_extremal_leves(neighbs, bias);
  }
  return neighbs;
}

}
