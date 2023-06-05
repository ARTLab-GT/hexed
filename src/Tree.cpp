#include <Tree.hpp>

namespace hexed
{

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

Tree* Tree::find_leaf(int ref_level, Eigen::VectorXi c, Eigen::VectorXi bias)
{
  return nullptr;
}

}
