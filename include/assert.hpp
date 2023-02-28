#ifndef HEXED_ASSERT_HPP_
#define HEXED_ASSERT_HPP_

#include <stdexcept>
#include <omp.h>
#include "config.hpp"

//! utilities for custom assertions
namespace hexed::assert
{

//! throws a `std::runtime_error` with message `message`, wrapped in a `#pragma omp critical` if necessary.
inline void throw_critical(const char* message)
{
  #if HEXED_THREADED
  if (omp_get_level()) {
    // if this is in a parallel region, only let one thread throw
    // or else instead of the error message all you see is "terminate called recursively"
    #pragma omp critical
    throw std::runtime_error(message);
  } else
  #endif
  {
    // if not in a parallel region, don't use the pragma cause that messes up try/catch
    throw std::runtime_error(message);
  }
}

}

/*! \brief Assert something with a helpful error message.
 * \details If `expression` is false, throws a `std::runtime_error` with a message
 * that includes `message` plus some additional info for debugging.
 * Works inside single-threaded regions as well as OpenMP parallel regions.
 */
#define HEXED_ASSERT(expression, message) \
  if (!(expression)) { \
    char buffer [1000]; \
    snprintf(buffer, 1000, "%s\n" \
                           "Exact cause: assertion `%s` failed in %s.\n" \
                           "Assertion invoked at line %d of %s in function %s.", \
             (message), #expression, __FUNCTION__, __LINE__, __FILE__, __PRETTY_FUNCTION__); \
    assert::throw_critical(buffer); \
  } \

#endif
