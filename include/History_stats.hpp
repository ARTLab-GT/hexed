#ifndef HEXED_HISTORY_STATS_HPP_
#define HEXED_HISTORY_STATS_HPP_

#include "math.hpp"

namespace hexed {

class History_stats {
  public:
  History_stats(double iteration_fraction = 0.1);
  void add_sample(Int iter, double value);
  inline Int n_sample() const {return _n_sample;}
  inline Int last_iter() const {return _last_iter;}
  inline double last_value() const {return _last_value;}
  inline double mean() const {return _mean;}
  inline double std_dev() const {return _std_dev;}
  inline double deriv() const {return _mean_deriv;}
  inline double deriv_std_dev() const {return _deriv_std_dev;}

  private:
  double _iter_frac;
  Int _last_iter;
  double _last_value;
  Int _n_sample;
  double _mean;
  double _variance;
  double _std_dev;
  double _mean_deriv;
  double _deriv_std_dev;
  Mat<2, 2> _lhs;
  Mat<2> _rhs;
};

}
#endif
