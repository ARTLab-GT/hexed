#ifndef HEXED_HISTORY_STATS_HPP_
#define HEXED_HISTORY_STATS_HPP_

#include "math.hpp"

namespace hexed {

/*! \brief Computes statistics about the convergence history of some variable.
 * \details Computes a quadratic curve fit to the entire convergence history
 * with a weighted least squares error metric.
 * The weights are exponentially decaying such that the last `iteration_fraction` iterations
 * should dominate the error metric.
 * This formulation only needs O(1) memory because each time a new sample point is added,
 * the new weighted normal equations can be computed directly from the previous normal equations
 * and the data from the latest sample point.
 */
class History_stats {
  public:
  //! \param iteration_fraction What fraction of the iterations should have a significant impact on the curve fit
  History_stats(double iteration_fraction);
  //! \brief Adds a new sample point and updates the curve fit to include it.
  //! \details `iter` must be nonnegative
  //! and strictly greater than the `iter` value passed to the last `add_sample` call.
  void add_sample(Int iter, double value);
  //! \brief Retrieves the `iteration_fraction` argument passed to the \ref History_stats "constructor".
  double iteration_fraction() const {return _iter_frac;}
  //! \brief The number of sample points added with `add_sample()` (Not the number of iterations!)
  inline Int n_sample() const {return _n_sample;}
  //! \brief The value of `iter` passed to the last `add_sample()` call
  inline Int last_iter() const {return _last_iters[0];}
  //! \brief The `value` passed to the last `add_sample()` call
  inline double last_value() const {return _last_values[0];}
  //! \brief The value of the curve fit at the last iteration (`last_iter()`).
  inline double smoothed() const {return _smoothed;}
  //! \brief The derivative of the curve fit with respect to iteration count at the last iteration (`last_iter()`)
  inline double trend() const {return _trend;}
  //! \brief The second derivative of the curve fit with respect to iteration count
  //! at the last iteration (`last_iter()`)
  inline double curvature() const {return _curve;}
  void reset(); //!< \brief Clears the sample history.

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
