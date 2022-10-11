#ifndef HEXED_ASSERT_HPP_
#define HEXED_ASSERT_HPP_

#include <stdexcept>
#include "config.hpp"

namespace hexed::assert
{

inline void throw_critical(const char* message)
{
  #pragma omp critical
  throw std::runtime_error(message);
}

}

#define HEXED_ASSERT(expression, message) \
  if (!(expression)) { \
    { \
      char buffer [1000]; \
      snprintf(buffer, 1000, "%s (in %s)", (message), __PRETTY_FUNCTION__); \
      assert::throw_critical(buffer); \
    } \
  } \

#endif
