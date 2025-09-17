#include <hexed/History_stats.hpp>
#include <hexed/Printer.hpp> //! \todo delete this

namespace hexed {

History_stats::History_stats(double iteration_fraction)
: _iter_frac{iteration_fraction}
, _last_iter{-1}
, _last_value{std::nan("")}
, _n_sample{0}
, _smoothed{std::nan("")}
, _trend{std::nan("")}
, _lhs{Mat<2, 2>::Zero()}
, _rhs{Mat<2>::Zero()}
{}

void History_stats::add_sample(Int iter, double value) {
  HEXED_ASSERT(iter >= 0, "Iteration values must be nonnegative.")
  HEXED_ASSERT(iter > _last_iter, "Iteration values must be increasing.")
  if (_n_sample > 0) {
    double new_weight = std::min(.5, (iter - _last_iter)/(_iter_frac*std::max<Int>(iter, 1)));
    double old_weight = math::pow(1 - new_weight, 2);
    new_weight *= new_weight;
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
      // shouldn't be aliasing because entries are traversed in matching order
      _lhs = old_weight*_lhs + new_weight*lhs_update;
      _rhs = old_weight*_rhs + new_weight*Mat<2>{iter*value, value};
    }
    Mat<2> soln = _lhs.inverse()*_rhs; // since this matrix is super small inverse is fine
    _smoothed = soln(0)*iter + soln(1);
    _trend = soln(0);
  } else {
    _smoothed = value;
  }
  _last_iter = iter;
  _last_value = value;
  ++_n_sample;
}

}
