#ifndef HEXED_CONVERGENCE_MONITOR_HPP_
#define HEXED_CONVERGENCE_MONITOR_HPP_

#include "History_stats.hpp"
#include "math.hpp"

namespace hexed {

class Convergence_monitor {
  public:
  Convergence_monitor(double iteration_fraction);
  void add_sample(Int iteration, double value);
  inline Int n_sample() const {return _value.n_sample();}
  inline Int last_iter() const {return _value.last_iter();}
  double iteration_fraction() const {return _value.iteration_fraction();}
  inline double last_value() const {return _value.last_value();}
  inline double smoothed() const {return _value.smoothed();}
  inline double trend() const {return _value.trend();}
  inline double curvature() const {return _value.curvature();}
  inline double noise() const {return std::max(0., _noise.smoothed());}
  inline double noise_trend() const {return _noise.trend();}
  inline double noise_curvature() const {return _noise.curvature();}
  bool converged(math::Tolerance trend_tol, math::Tolerance noise_tol) const;

  private:
  History_stats _value;
  History_stats _noise;
};

}
#endif
