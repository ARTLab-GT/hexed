#ifndef CARTDG_STOPWATCH_HPP_
#define CARTDG_STOPWATCH_HPP_

#include <chrono>

namespace cartdg
{

/*
 * A class for conveniently measuring execution time. Can be started
 * and stopped multiple times. Thows an exception if started twice without
 * stopping, or vice-versa. `time` returns the cumulative elapsed time
 * between `start`s and `pause`s. `n_calls` returns the number of `start`/`pause`
 * cycles performed. `reset` resets `time` and `n_calls` to 0.
 */
class Stopwatch
{
  int n = 0;
  double t = 0.;
  bool r = false;
  std::chrono::high_resolution_clock::time_point time_started;

  public:
  void start();
  void pause();
  void reset();
  bool running();
  int n_calls();
  double time();
};

}
#endif
