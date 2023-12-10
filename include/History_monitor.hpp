#ifndef HEXED_HISTORY_MONITOR_HPP_
#define HEXED_HISTORY_MONITOR_HPP_

#include <vector>

namespace hexed
{

/*! \brief monitors the history of some variable over iterations and computes the bounds over a specified window
 * \details The class user supplies the value of the desired variable at certain iterations.
 * The `History_monitor` will compute the maximum and minimum of these values over the last some-percent of the iterations.
 * In case the number of iterations becomes large, there is a user-specified maximum number of samples,
 * and the `History_monitor` will only record a fraction of the supplied data points at a frequency chosen to maintain the specified buffer size.
 * Of course, values outside of the window are also forgotten.
 */
class History_monitor
{
  std::vector<int> _iterations;
  std::vector<double> _values;
  int _sz;
  int _start;
  int _end;

  public:
  /*!\param window_size Window size as a fraction of the iteration count.
   *   E.g. if `window_size = .3`, bounds will be computed over the last 30% of the iterations.
   * \param max_samples Maximum sample buffer size.
   */
  History_monitor(double window_size, int max_samples);
  void add_sample(int iteration, double value); //!< \brief stipulates that the value of the variable to be monitored is `value` at iteration `iteration`
  double max(); //!< \brief obtains the maximum of the variable over the window
  double min(); //!< \brief obtains the minimum of the variable over the window
};

}
#endif
