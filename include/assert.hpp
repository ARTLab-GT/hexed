#ifndef HEXED_ASSERT_HPP_
#define HEXED_ASSERT_HPP_

#include <stdexcept>
#include <omp.h>
#include "config.hpp"

//! \file assert.hpp utilities for custom assertions

//! utilities for custom assertions
namespace hexed::assert
{

class Exception : public std::exception
{
  std::string msg;
  public:
  inline Exception(std::string message) : msg{message} {}
  inline const char* what() const noexcept override {return msg.c_str();}
};

//! represents a fatal problem in the numerics of the code (such as nonphysical values)
//! as opposed to, for example, an out-of-bounds error or user error
//! \see \ref numerical_errors
class Numerical_exception : public Exception
{
  public:
  Numerical_exception(std::string message) : Exception(message) {}
};

//! throws a `std::runtime_error` with message `message`, wrapped in a `#pragma omp critical` if necessary.
//! Used in \ref HEXED_ASSERT
template <typename except_t = std::runtime_error>
void throw_critical(const char* message)
{
  #if HEXED_THREADED
  if (omp_get_level()) {
    // if this is in a parallel region, only let one thread throw
    // or else instead of the error message all you see is "terminate called recursively"
    #pragma omp critical
    throw except_t(message);
  } else
  #endif
  {
    // if not in a parallel region, don't use the pragma cause that messes up try/catch
    throw except_t(message);
  }
}

}

/*! \brief Assert something with a helpful error message.
 * \details If `expression` is false, throws an exception with a message
 * that includes `message` plus some additional info for debugging.
 * Works inside single-threaded regions as well as OpenMP parallel regions.
 * If desired, supply the type of exception as the third argument.
 * Exception type must be constructible from a string.
 * Defaults to `std::runtime_error`.
 */
#define HEXED_ASSERT(expression, message, ...) { \
  if (!(expression)) { \
    char buffer [1000]; \
    snprintf(buffer, 1000, "%s\n" \
                           "technical details: assertion `%s` failed in `%s`.\n" \
                           "Assertion invoked at line %d of %s in function %s.", \
             std::string(message).c_str(), #expression, __FUNCTION__, __LINE__, __FILE__, __PRETTY_FUNCTION__); \
    assert::throw_critical<__VA_ARGS__>(buffer); \
  } \
}

#endif
