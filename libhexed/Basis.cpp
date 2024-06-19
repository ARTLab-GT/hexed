#include <hexed/Basis.hpp>

namespace hexed
{

double Basis::max_cfl() const
{
  return -2*quadratic_safety()/min_eig_convection();
}

double Basis::step_ratio() const
{
  return .5/quadratic_safety();
}

Basis::Basis(int row_size_arg) : row_size(row_size_arg) {}
Basis::~Basis() {}

Mat<> Basis::nodes() const
{
  Mat<> n(row_size);
  for (int i_node = 0; i_node < row_size; ++i_node) n(i_node) = node(i_node);
  return n;
}

Mat<dyn, dyn> Basis::interpolate(const Mat<>& sample) const
{
  Mat<dyn, dyn> interp {Mat<dyn, dyn>::Ones(sample.size(), row_size)};
  for (int i_node = 0; i_node < row_size; ++i_node) {
    for (int j_node = 0; j_node < row_size; ++j_node) {
      if (i_node != j_node) {
        interp.col(i_node).array() *= (sample - Mat<>::Constant(sample.size(), node(j_node))).array()
                                      /(node(i_node) - node(j_node));
      }
    }
  }
  return interp;
}

Mat<dyn, dyn> Basis::restrict(int i_half) const
{
  return interpolate(.5*(nodes() + Mat<>::Constant(row_size, double(i_half))));
}

}
