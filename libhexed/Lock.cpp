#include <Lock.hpp>

namespace hexed {

Lock::Set::Set(Lock& ref)
: _lock{&ref}
{
  #if HEXED_THREADED
  omp_set_nest_lock(&_lock->_l);
  #endif
}

Lock::Set::Set(Lock::Set&& that)
: _lock{nullptr}
{
  *this = std::move(that);
}

Lock::Set& Lock::Set::operator=(Lock::Set&& that) {
  #if HEXED_THREADED
  if (_lock) omp_unset_nest_lock(&_lock->_l);
  #endif
  _lock = that._lock;
  that._lock = nullptr;
  return *this;
}


Lock::Set::~Set() {
  #if HEXED_THREADED
  if (_lock) omp_unset_nest_lock(&_lock->_l);
  #endif
}

Lock::Lock() {
  #if HEXED_THREADED
  omp_init_nest_lock(&_l);
  #endif
}

Lock::Lock(const Lock&) : Lock() {}
void Lock::operator=(const Lock&) {}

Lock::~Lock() {
  #if HEXED_THREADED
  omp_destroy_nest_lock(&_l);
  #endif
}

std::optional<Lock::Set> Lock::test() {
  std::optional<Set> set;
  #if HEXED_THREADED
  if (omp_test_nest_lock(&_l)) { // Increments nesting count. We now own the lock, so we can set it without blocking.
    set.emplace(*this); // increments nesting count again
    omp_unset_nest_lock(&_l); // Decrement nesting count to undo initial test. We still own the lock.
  }
  #endif
  return set;
}

}
