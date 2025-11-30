#include <hexed/History_stats.hpp>

namespace hexed {

History_stats::History_stats(double iter_frac)
: _iter_frac{iter_frac}
, _last_iters{-1, -1}
, _last_values{std::nan(""), std::nan("")}
, _n_sample{0}
, _smoothed{std::nan("")}
, _trend{std::nan("")}
, _curve{std::nan("")}
, _lhs{Mat<3, 3>::Zero()}
, _rhs{Mat<3>::Zero()}
{}

void History_stats::add_sample(Int iter, double value) {
  HEXED_ASSERT(iter >= 0, "Iteration values must be nonnegative.")
  HEXED_ASSERT(iter > _last_iters[0], "Iteration values must be increasing.")
  if (!std::isfinite(value)) return;
  if (_n_sample > 1) {
    double new_weight = std::min(.5, (iter - _last_iters[0])/(_iter_frac*std::max<Int>(iter, 1)));
    double old_weight = math::pow(1 - new_weight, 2);
    new_weight *= new_weight;
    if (_n_sample == 2) {
      // matrix and rhs to fit the last 3 points exactly
      _lhs << _last_iters[1]*_last_iters[1]*.5, _last_iters[1], 1.,
              _last_iters[0]*_last_iters[0]*.5, _last_iters[0], 1.,
                                  iter*iter*.5,           iter, 1.;
      _rhs << _last_values[1], _last_values[0], value;
      // change them into the normal equations
      _rhs = _lhs.transpose()*_rhs;
      _lhs = _lhs.transpose()*_lhs; // shouldn't be aliasing because matrix multiplication automatically copies
    } else {
      // compute the matrix for the normal equations for the last point
      Mat<1, 3> last_point {.5*iter*iter, double(iter), 1.};
      Mat<3, 3> lhs_update = last_point.transpose()*last_point;
      // use this matrix to update the global normal equations.
      // there shouldn't be aliasing because entries are traversed in matching order
      _lhs = old_weight*_lhs + new_weight*lhs_update;
      _rhs = old_weight*_rhs + new_weight*last_point.transpose()*value;
    }
    Mat<3> soln = _lhs.inverse()*_rhs; // since this matrix is super small inverse is fine
    _smoothed = .5*soln(0)*iter*iter + soln(1)*iter + soln(2);
    _trend = soln(0)*iter + soln(1);
    _curve = soln(0);
  } else {
    _smoothed = value;
  }
  _last_iters[1] = _last_iters[0];
  _last_iters[0] = iter;
  _last_values[1] = _last_values[0];
  _last_values[0] = value;
  ++_n_sample;
}

}
