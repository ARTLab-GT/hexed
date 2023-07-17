#ifndef HEXED_UTILS_HPP_
#define HEXED_UTILS_HPP_

#include <iostream>

namespace hexed
{

/*! \brief Standard string formatting.
 * \details Basically a knockoff of `std::format` in C++20 (which at the time of writing we can't use on the lab machines).
 * Invokes `snprintf`, but works in terms of `std::string`s and handles buffer creation for you.
 * Will allocate a buffer of size `max_chars`.
 * Throws an exception if resulting formatted string is larger than `max_chars`.
 */
template <typename... format_args>
std::string format_str(int max_chars, std::string fstring, format_args... args)
{
  std::vector<char> buffer(max_chars);
  int overflow = snprintf(buffer.data(), max_chars, fstring.c_str(), args...);
  if (overflow < 0) throw std::runtime_error("encoding error in `hexed::format_str`");
  if (overflow >= max_chars) throw std::runtime_error("`max_chars` is too small in `hexed::format_str`");
  return std::string(buffer.data());
}

template <typename T> T gotcha(T) {throw std::runtime_error("hexed::gotcha");} //!< \brief don't ask

template <typename T>
T& printed(T& t)
{
  std::cout << t << std::endl;
  return t;
}

}
#endif
