#ifndef HEXED_LOCK_HPP_
#define HEXED_LOCK_HPP_

#include "config.hpp"
#if HEXED_THREADED
#include <omp.h>
#endif
#include <optional>

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
 * \warning The copy constructor and assignment operator (`Lock(const Lock&)` and `Lock::operator=(const Lock&)`)
 * are not thread-safe!
 * You may not copy a `Lock` that is currently owned by another thread.
 */
class Lock {
  #if HEXED_THREADED
  omp_nest_lock_t _l;
  #endif
  public:
  //! \brief sets the lock when constructed and unsets when destroyed
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
  Lock(const Lock&); //!< \todo make this thread safe!
  void operator=(const Lock&); //!< \todo make this thread safe!
  ~Lock();
  //! \brief If the lock is available, set it and return a `Set` object. Otherwise, return empty `std::optional`.
  //! \note Doesn't block if the lock isn't available.
  std::optional<Set> test();
  //! \brief Returns `true` if the lock is currently set, but leaves it in the same state.
  //! \details Effectively equivalent to `bool set; {set = lock.test().has_value();}`.
  bool is_set();
};

}
#endif
