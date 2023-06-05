#include <Tree.hpp>

namespace hexed
{

Tree::Tree(int nd, double root_size, Mat<> origin)
: root_sz{root_size}, ref_level{0}, coords{Eigen::VectorXi::Zero(nd)}, n_dim{nd}
{
  HEXED_ASSERT(origin.size() >= n_dim, "`origin` is too small");
  orig = origin(Eigen::seqN(0, n_dim));
}

Mat<> Tree::origin() const {return orig;}
int Tree::refinement_level() const {return ref_level;}
Eigen::VectorXi Tree::coordinates() const {return coords;}
double Tree::nominal_size() const {return root_sz;}
Mat<> Tree::nominal_position() const {return coords.cast<double>();} // note: wrong

}
