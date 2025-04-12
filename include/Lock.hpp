#ifndef HEXED_LOCK_HPP_
#define HEXED_LOCK_HPP_

#include <omp.h>
#include "config.hpp"

namespace hexed {

/*! \brief wrapper for [OpenMP **nested** lock routines](https://www.openmp.org/spec-html/5.0/openmpse31.html).
 * \details This class can be used to prevent data races by OpenMP threads.
 * It contains an OpenMP lock variable.
 * To set the lock, construct a `Set` object from it.
 * When the `Set` object is destroyed, the lock will be unset.
 * If Hexed is not compiled with OpenMP, this class does nothing.
 * Use like this:
 * ~~~
 * Lock l;
 * #pragma omp parallel for
 * for (int i = 0; i < N; ++i) {
 *   Lock::Set s(l); // sets lock
 *   // only one thread at a time can execute any statements here
 * } // lock is released because `s` is destroyed
 * ~~~
 */
class Lock {
  #if HEXED_THREADED
  omp_nest_lock_t _l;
  #endif
  public:
  //! sets the lock when constructed and unsets when destroyed
  class Set {
    public:
    Set(Lock&);
    Set(Set&& that);
    Set(const Set&) = delete;
    Set& operator=(Set&& that);
    Set& operator=(const Set&) = delete;
    ~Set();
    private:
    Lock* _lock;
  };
  Lock();
  ~Lock();
};

}
#endif
