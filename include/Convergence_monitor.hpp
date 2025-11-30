#ifndef HEXED_CONVERGENCE_MONITOR_HPP_
#define HEXED_CONVERGENCE_MONITOR_HPP_

#include "History_stats.hpp"
#include "math.hpp"

namespace hexed {

/*! \brief Monitors whether some variable has converged to a definite value.
 * \details This is essentially a package of two `History_stats` objects,
 * one of which monitors the history of the variable itself,
 * and the other of which monitors the history of the deviation of the variable from the curve fit,
 * which is referred to as the _noise_.
 * The variable is considered to be converged if there is no significant drift over time
 * (the `trend()` and `curvature()` are sufficiently small)
 * and _either_ of the following is true:
 * - The `noise()` is sufficiently small, indicating that the variable has converged to a single value.
 * - The `noise()` has stagnated (its `trend()` and `curvature()` are sufficiently small),
 *   indicating that the variable is displaying an oscillatory pattern that is not changing over time,
 *   and running more iterations will not provide any new information.
 */
class Convergence_monitor {
  public:
  //! \param iteration_fraction The `iteration_fraction` passed to the \ref History_stats::History_stats "History_stats"
  Convergence_monitor(double iteration_fraction);
  //! \brief Adds a sample to the value and noise monitors.
  //! \see `History_stats::add_sample`
  void add_sample(Int iteration, double value);
  //! \brief Forwards to `History_stats::n_sample`
  inline Int n_sample() const {return _value.n_sample();}
  //! \brief Forwards to `History_stats::last_iter`
  inline Int last_iter() const {return _value.last_iter();}
  //! \brief Forwards to `History_stats::iteration_fraction`
  double iteration_fraction() const {return _value.iteration_fraction();}
  //! \brief Forwards to `History_stats::last_value` for the value `History_stats`
  inline double last_value() const {return _value.last_value();}
  //! \brief Forwards to `History_stats::smoothed` for the value `History_stats`
  inline double smoothed() const {return _value.smoothed();}
  //! \brief Forwards to `History_stats::trend` for the value `History_stats`
  inline double trend() const {return _value.trend();}
  //! \brief Forwards to `History_stats::curvature` for the value `History_stats`
  inline double curvature() const {return _value.curvature();}
  //! \brief Forwards to `History_stats::smoothed` for the noise `History_stats`
  inline double noise() const {return std::max(0., _noise.smoothed());}
  //! \brief Forwards to `History_stats::trend` for the noise `History_stats`
  inline double noise_trend() const {return _noise.trend();}
  //! \brief Forwards to `History_stats::curvature` for the noise `History_stats`
  inline double noise_curvature() const {return _noise.curvature();}
  /*! \brief Checks whether the variable has converged to within the specified tolerances.
   * \details `trend_tol` applies to the variation of the curve fit over the last `iteration_fraction` iterations,
   * where the absolute value of the linear and quadratic components are taken separately
   * so that they can't cancel out.
   * `trend_tol.rel` is relative to `smoothed() + noise()` so that it can still be applicable to quantities
   * that converge to zero if they do so in an oscillatory manner.
   * `noise_tol.rel` is relative to `smoothed()`,
   * and `trend_tol` is used to assess whether the noise has stagnated,
   * with `trend_tol.rel` being relative to `noise()` in this case.
   */
  bool converged(math::Tolerance trend_tol, math::Tolerance noise_tol) const;

  private:
  History_stats _value;
  History_stats _noise;
};

}
#endif
