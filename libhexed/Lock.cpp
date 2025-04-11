#include <Lock.hpp>

namespace hexed {

Lock::Acquire::Acquire(Lock& ref)
: _lock{&ref}
{
  #if HEXED_THREADED
  omp_set_nest_lock(&_lock->_l);
  #endif
}

Lock::Acquire::Acquire(Lock::Acquire&& that)
: _lock{nullptr}
{
  *this = std::move(that);
}

Lock::Acquire& Lock::Acquire::operator=(Lock::Acquire&& that) {
  #if HEXED_THREADED
  if (_lock) omp_unset_nest_lock(&_lock->_l);
  #endif
  _lock = that._lock;
  that._lock = nullptr;
  return *this;
}


Lock::Acquire::~Acquire() {
  #if HEXED_THREADED
  if (_lock) omp_unset_nest_lock(&_lock->_l);
  #endif
}

Lock::Lock() {
  #if HEXED_THREADED
  omp_init_nest_lock(&_l);
  #endif
}

Lock::~Lock() {
  #if HEXED_THREADED
  omp_destroy_nest_lock(&_l);
  #endif
}

}
