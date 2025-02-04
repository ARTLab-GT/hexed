#include <Lock.hpp>

namespace hexed {

Lock::Acquire::Acquire(Lock& ref)
: lock{ref}
{
  #if HEXED_THREADED
  omp_set_nest_lock(&lock.l);
  #endif
}

Lock::Acquire::~Acquire() {
  #if HEXED_THREADED
  omp_unset_nest_lock(&lock.l);
  #endif
}

Lock::Lock() {
  #if HEXED_THREADED
  omp_init_nest_lock(&l);
  #endif
}

Lock::~Lock() {
  #if HEXED_THREADED
  omp_destroy_nest_lock(&l);
  #endif
}

}
