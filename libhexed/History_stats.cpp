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
  UPDATE(_variance, math::pow(value - _mean, 2), 1)
  UPDATE(_mean_deriv, (_mean - old_mean)/(iter - _last_iter), 1)
  double old_std_dev = _std_dev;
  _std_dev = std::sqrt(_variance);
  #undef UPDATE
  _last_iter = iter;
  _last_value = value;
  ++_n_sample;
}

}
