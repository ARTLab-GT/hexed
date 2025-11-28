#ifndef HEXED_HISTORY_STATS_HPP_
#define HEXED_HISTORY_STATS_HPP_

#include "math.hpp"

namespace hexed {

class History_stats {
  public:
  History_stats(double iteration_fraction);
  void add_sample(Int iter, double value);
  double iteration_fraction() const {return _iter_frac;}
  inline Int n_sample() const {return _n_sample;}
  inline Int last_iter() const {return _last_iters[0];}
  inline double last_value() const {return _last_values[0];}
  inline double smoothed() const {return _smoothed;}
  inline double trend() const {return _trend;}
  inline double curvature() const {return _curve;}

  private:
  double _iter_frac;
  Int _last_iters [2];
  double _last_values [2];
  Int _n_sample;
  double _smoothed;
  double _trend;
  double _curve;
  Mat<3, 3> _lhs;
  Mat<3> _rhs;
};

}
#endif
