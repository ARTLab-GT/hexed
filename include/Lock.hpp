#ifndef HEXED_LOCK_HPP_
#define HEXED_LOCK_HPP_

#include <omp.h>
#include "config.hpp"

namespace hexed
{

/*! \brief wrapper for [OpenMP lock routines](https://www.openmp.org/spec-html/5.0/openmpse31.html).
 * \details This class can be used to prevent data races by OpenMP threads.
 * It contains an OpenMP lock variable.
 * To acquire (aka set) the lock, construct an `Acquire` object from it.
 * When the `Acquire` object is destroyed, the lock will be released (unset).
 * If Hexed is not compiled with OpenMP, this class does nothing.
 * Use like this:
 * ~~~
 * Lock l;
 * #pragma omp parallel for
 * for (int i = 0; i < N; ++i) {
 *   Lock::Acquire a(l); // acquires lock
 *   // only one thread at a time can execute any statements here
 * } // lock is released because `a` is destroyed
 * ~~~
 */
class Lock
{
  #if HEXED_THREADED
  omp_lock_t l;
  #endif
  public:
  //! acquires the lock when constructed and releases when destroyed
  class Acquire
  {
    Lock& lock;
    public:
    Acquire(Lock&);
    ~Acquire();
  };
  Lock();
  ~Lock();
};

}
#endif
