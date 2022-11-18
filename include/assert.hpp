#ifndef HEXED_ASSERT_HPP_
#define HEXED_ASSERT_HPP_

#include <stdexcept>
#include <omp.h>
#include "config.hpp"

namespace hexed::assert
{

inline void throw_critical(const char* message)
{
  if (omp_get_level()) {
    // if this is in a parallel region, only let one thread throw
    // or else instead of the error message all you see is "terminate called recursively"
    #pragma omp critical
    throw std::runtime_error(message);
  } else {
    // if not in a parallel region, don't use the pragma cause that messes up try/catch
    throw std::runtime_error(message);
  }
}

}

#define HEXED_ASSERT(expression, message) \
  if (!(expression)) { \
    char buffer [1000]; \
    snprintf(buffer, 1000, "%s (in %s)", (message), __PRETTY_FUNCTION__); \
    assert::throw_critical(buffer); \
  } \

#endif
