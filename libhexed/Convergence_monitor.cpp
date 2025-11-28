#include <hexed/Convergence_monitor.hpp>

namespace hexed {

Convergence_monitor::Convergence_monitor(double iter_frac)
: _value(iter_frac)
, _noise(iter_frac)
{}

void Convergence_monitor::add_sample(Int iter, double value) {
  _value.add_sample(iter, value);
  _noise.add_sample(iter, 2*std::abs(value - _value.smoothed()));
}

bool Convergence_monitor::converged(math::Tolerance trend_tol, math::Tolerance noise_tol) const {
  double trend_tol_total = trend_tol.rel*(std::abs(smoothed()) + noise()) + trend_tol.abs;
  double iter_factor = iteration_fraction()*last_iter();
  bool trend_conv = std::abs(trend()*iter_factor) + std::abs(.5*curvature()*iter_factor*iter_factor) < trend_tol_total;
  double noise_tol_total = noise_tol.rel*std::abs(smoothed()) + noise_tol.abs;
  bool noise_conv = noise() + .1*std::abs(noise_trend()) < noise_tol_total;
  bool noise_trend_conv = std::abs(noise_trend())*iter_factor + std::abs(.5*noise_curvature()*iter_factor*iter_factor)
                          < trend_tol.rel*noise() + trend_tol.abs;
  return trend_conv && (noise_conv || noise_trend_conv);
}

}
