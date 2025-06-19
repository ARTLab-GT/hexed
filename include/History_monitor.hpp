#ifndef HEXED_HISTORY_MONITOR_HPP_
#define HEXED_HISTORY_MONITOR_HPP_

#include <vector>
#include "math.hpp"

namespace hexed {

/*! \brief monitors the history of some variable over iterations and computes the bounds over a specified window
 * \details The class user supplies the value of the desired variable at certain iterations.
 * The `History_monitor` will compute the maximum and minimum of these values
 * over the last some-percent of the iterations.
 * In case the number of iterations becomes large, there is a user-specified maximum number of samples,
 * and the `History_monitor` will only record a fraction of the supplied data points
 * at a frequency chosen to maintain the specified buffer size.
 * Of course, values outside of the window are also forgotten.
 */
class History_monitor {
  Int _samples; // size of sample storage vector
  Int _min_samples; // minimum number of samples required for conclusive min/max
  Int _start; // index of the oldest sample
  Int _sz; // current number of samples stored
  // Samples will be stored in the following vectors.
  // Since the max number of samples is known, to allow efficient addition/removal of items at the ends of the vectors,
  // the vectors are allocated with `_samples` elements.
  // Not all the elements may be filled with samples at any given time and periodic indexing is used.
  // At any given time, the current set of samples goes from index `_start` to `_start + _sz`, wrapping around
  // from the end to the beginning if necessary.
  std::vector<int> _iterations; // vector of the iterations of all the samples stored
  std::vector<double> _values; // vector of the values of all the samples stored
  double _win_sz; // fraction of the total iterations to store samples for
  double _add_threshold; // the next sample stored will be the first one whose iteration exceeds this value
  // current values of min and max of recorded samples
  double _min;
  double _max;

  public:
  /*!\param window_size Window size as a fraction of the iteration count.
   *   E.g. if `window_size = .3`, bounds will be computed over the last 30% of the iterations.
   * \param max_samples Maximum sample buffer size.
   * \param min_samples `min()` and `max()` will give return +- large numbers until this many samples have been stored
   *   (and only samples within the `window_size` will be stored).
   */
  History_monitor(double window_size, Int max_samples, Int min_samples = 2);
  //! \brief stipulates that the value of the variable to be monitored is `value` at iteration `iteration`
  void add_sample(Int iteration, double value);
  double min() const; //!< \brief obtains the minimum of the variable over the window
  double max() const; //!< \brief obtains the maximum of the variable over the window
  void clear(); //!< \brief deletes the entire history of the monitor as though it had been constructed again
};

}
#endif
