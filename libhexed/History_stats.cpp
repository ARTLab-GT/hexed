#include <hexed/History_stats.hpp>
#include <hexed/Printer.hpp> //! \todo delete this

namespace hexed {

History_stats::History_stats(double iteration_fraction)
: _iter_frac{iteration_fraction}
, _last_iter{-1}
, _last_value{std::nan("")}
, _n_sample{0}
, _mean{std::nan("")}
, _variance{std::nan("")}
, _std_dev{std::nan("")}
, _mean_deriv{std::nan("")}
, _deriv_std_dev{std::nan("")}
, _lhs{Mat<2, 2>::Zero()}
, _rhs{Mat<2>::Zero()}
{}

void History_stats::add_sample(Int iter, double value) {
  HEXED_ASSERT(iter >= 0, "Iteration values must be nonnegative.")
  HEXED_ASSERT(iter > _last_iter, "Iteration values must be increasing.")
  #define UPDATE(field, expr, sample) { \
    if (_n_sample == sample) { \
      field = expr; \
    } \
    if (_n_sample > sample) { \
      field += ((expr) - field)*std::min(1., (iter - _last_iter)/(_iter_frac*iter)); \
    } \
  }
  double old_mean = _mean;
  UPDATE(_mean, value, 0)
  UPDATE(_mean_deriv, (_mean - old_mean)/(iter - _last_iter), 1)
  if (_n_sample > 0) {
    double new_weight = std::min(_iter_frac, (iter - _last_iter)/(_iter_frac*iter));
    if (_n_sample == 1) {
      _lhs << _last_iter, 1.,
                    iter, 1.;
      _rhs << _last_value, value;
      _rhs = _lhs.transpose()*_rhs;
      _lhs = _lhs.transpose()*_lhs; // shouldn't be aliasing because matrix multiplication automatically copies
    } else {
      Mat<2, 2> lhs_update;
      lhs_update << iter*iter, iter,
                         iter,   1.;
      _lhs += new_weight*(lhs_update - _lhs); // shouldn't be aliasing because entries are traversed in matching order
      _rhs += new_weight*(Mat<2>{iter*value, value} - _rhs);
    }
    Mat<2> soln = _lhs.inverse()*_rhs; // since this matrix is super small inverse is fine
    _mean = soln(0)*iter + soln(1);
    _mean_deriv = soln(0);
  }
  UPDATE(_variance, math::pow(value - _mean, 2), 1)
  double old_std_dev = _std_dev;
  _std_dev = std::sqrt(_variance);
  UPDATE(_deriv_std_dev, (_std_dev - old_std_dev)/(iter - _last_iter), 2)
  #undef UPDATE
  _last_iter = iter;
  _last_value = value;
  ++_n_sample;
}

}
