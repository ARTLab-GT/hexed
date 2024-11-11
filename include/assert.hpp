#ifndef HEXED_ASSERT_HPP_
#define HEXED_ASSERT_HPP_

#include <stdexcept>
#include <vector>
#include <string>
#include <omp.h>
#include "config.hpp"

//! \file assert.hpp utilities for custom assertions

namespace hexed {

/*! \brief Standard string formatting.
 * \details Basically a knockoff of `std::format` in C++20 (which at the time of writing we can't use on the lab machines).
 * Invokes `snprintf`, but works in terms of `std::string`s and handles buffer creation for you.
 * Will allocate a buffer of size `max_chars`.
 * Throws an exception if resulting formatted string is larger than `max_chars`.
 * \headerfile utils.hpp
 */
template <typename... format_args>
std::string format_str(int max_chars, std::string fstring, format_args... args) {
  std::vector<char> buffer(max_chars);
  int overflow = snprintf(buffer.data(), max_chars, fstring.c_str(), args...);
  if (overflow < 0) throw std::runtime_error("encoding error in `hexed::format_str`");
  if (overflow >= max_chars) throw std::runtime_error("`max_chars` is too small in `hexed::format_str`");
  return std::string(buffer.data());
}

//! \brief utilities for custom assertions
namespace assert {

class Exception : public std::exception {
  public:
  inline Exception(std::string m) : _msg{m} {}
  virtual std::string name() const = 0;
  std::string message() const {return _msg;}
  inline const char* what() const noexcept override {return _msg.c_str();}
  private:
  std::string _msg;
};

//! \brief Represents a fatal problem in the numerics of the code (such as nonphysical values)
//! \details as opposed to, for example, an out-of-bounds error or user error
//! \see \ref numerical_error
class Numerical_exception : public Exception {
  public:
  inline std::string name() const override {return "Numerical exception";}
  Numerical_exception(std::string message) : Exception(message) {}
};

//! \brief Represents an exception which clearly results from a mistake made by the user.
class User_error : public Exception {
  public:
  inline std::string name() const override {return "User error";}
  inline User_error(std::string message) : Exception(message) {}
};

//! \brief Indicates that the user invoked functionality which should be implemented in the future but isn't yet.
class Not_implemented_error : public Exception {
  public:
  inline std::string name() const override {return "Feature-not-implemented error";}
  inline Not_implemented_error(std::string message) : Exception(message) {}
};

//! \brief Indicates that a situation has been encountered which should not be possible, regardless of user input.
//! \details In other words, a bug has been caught.
class Internal_error : public Exception {
  public:
  inline std::string name() const override {return "Internal error";}
  inline Internal_error(std::string message) : Exception(message) {}
};

//! throws a `std::runtime_error` with message `message`, wrapped in a `#pragma omp critical` if necessary.
//! Used in \ref HEXED_ASSERT
template <typename except_t = Internal_error>
void throw_critical(const char* message) {
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

} // namespace assert
} // namespace hexed

/*! \brief Throws an exception with an informative error message.
 * \details Throws an exception with a message
 * that includes `message` plus some additional info for debugging.
 * Works inside single-threaded regions as well as OpenMP parallel regions.
 * If desired, supply the type of exception as the third argument.
 * Exception type must be constructible from a string.
 * Defaults to `std::runtime_error`.
 */
#define HEXED_THROW(message, ...) { \
  hexed::assert::throw_critical<__VA_ARGS__>(hexed::format_str(1000, \
    "%s\n" \
    "  At: line %d of `%s`\n" \
    "  In: %s", \
    std::string(message).c_str(), __LINE__, __FILE__, __PRETTY_FUNCTION__).c_str()); \
}

/*! \brief Assert something with an informative error message.
 * \details If `expression` is false, throws an exception with `HEXED_THROW`.
 * `message` and an optional third argument are passed to `HEXED_THROW`.
 */
#define HEXED_ASSERT(expression, message, ...) { \
  if (!(expression)) { \
    HEXED_THROW(hexed::format_str(1000, \
      "%s\n" \
      "  Assertion `%s` failed in `%s`.\n", \
      std::string(message).c_str(), #expression, __FUNCTION__) \
      __VA_OPT__(,) __VA_ARGS__ \
    ); \
  } \
}

#endif
