#ifndef HEXED_STOPWATCH_HPP_
#define HEXED_STOPWATCH_HPP_

#include <chrono>

namespace hexed
{

/*! \brief A class for conveniently measuring execution time.
 * \details Measures wall clock time.
 * Can be started and stopped multiple times.
 */
class Stopwatch
{
  int n = 0;
  double t = 0.;
  bool r = false;
  std::chrono::high_resolution_clock::time_point time_started;

  public:
  //! RAII-style operation of a stopwatch
  class Operator
  {
    Stopwatch& sw;
    public:
    Operator(Stopwatch& stopwatch) : sw{stopwatch} {sw.start();} //!< starts the stopwatch
    ~Operator() {sw.pause();} //!< stops the stopwatch
  };

  //! \brief Starts measuring time.
  //! \details Throws an exception if already running.
  void start();
  //! \brief pauses measurement and updates `time()`
  //! \details Throws an exception if not running.
  void pause();
  void reset(); //!< resets `time()` to zero
  bool running() const; //!< `false` if the stopwatch has been `stop`ped since the last time it was `start`ed
  int n_calls() const; //!< number of times the stopwatch has been `stop`ped.
  double time() const; //!< total elapsed time between `start()`s and `stop()`s since the last `reset()`
};

}
#endif
