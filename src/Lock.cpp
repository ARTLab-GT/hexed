#include <Lock.hpp>

namespace hexed
{

Lock::Acquire::Acquire(Lock& ref)
: lock{ref}
{
  omp_set_lock(&lock.l);
}

Lock::Acquire::~Acquire()
{
  omp_unset_lock(&lock.l);
}

Lock::Lock()
{
  omp_init_lock(&l);
}

Lock::~Lock()
{
  omp_destroy_lock(&l);
}

}
