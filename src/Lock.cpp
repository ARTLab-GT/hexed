#include <Lock.hpp>

namespace hexed
{

Lock::Acquire::Acquire(Lock& ref)
: lock{ref}
{
  #if HEXED_THREADED
  omp_set_lock(&lock.l);
  #endif
}

Lock::Acquire::~Acquire()
{
  #if HEXED_THREADED
  omp_unset_lock(&lock.l);
  #endif
}

Lock::Lock()
{
  #if HEXED_THREADED
  omp_init_lock(&l);
  #endif
}

Lock::~Lock()
{
  #if HEXED_THREADED
  omp_destroy_lock(&l);
  #endif
}

}
