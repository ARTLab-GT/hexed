#ifndef HEXED_STOPWATCH_HPP_
#define HEXED_STOPWATCH_HPP_

#include <chrono>

namespace hexed {

/*! \brief A class for conveniently measuring execution time.
 * \details Measures wall clock time.
 * Can be started and stopped multiple times.
 */
class Stopwatch {
  int n = 0;
  double t = 0.;
  bool r = false;
  std::chrono::steady_clock::time_point time_started;

  public:
  //! \brief RAII-style operation of a stopwatch
  class Operator {
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
  void reset(); //!< \brief resets `time()` to zero
  bool running() const; //!< \brief `false` if the stopwatch has been `stop`ped since the last time it was `start`ed
  int n_calls() const; //!< \brief number of times the stopwatch has been `stop`ped.
  double time() const; //!< \brief total elapsed time between `start()`s and `stop()`s since the last `reset()`
  Stopwatch operator+(Stopwatch other) const; //! adds `time` and `n_calls`, throws if either is running
  Stopwatch& operator+=(Stopwatch other);
};

}
#endif
