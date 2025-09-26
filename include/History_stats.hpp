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
  inline Int last_iter() const {return _last_iter;}
  inline double last_value() const {return _last_value;}
  inline double smoothed() const {return _smoothed;}
  inline double trend() const {return _trend;}

  private:
  double _iter_frac;
  Int _last_iter;
  double _last_value;
  Int _n_sample;
  double _smoothed;
  double _trend;
  Mat<2, 2> _lhs;
  Mat<2> _rhs;
};

}
#endif
