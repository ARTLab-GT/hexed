#ifndef HEXED_UTILS_HPP_
#define HEXED_UTILS_HPP_

#include <iostream>
#include <vector>

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

//! \brief prints its argument and then returns it
//! \details you can wrap this around an expression to print it without computing it again
template <typename T>
T& printed(T& t)
{
  std::cout << t << std::endl;
  return t;
}

template <typename T> T* new_copy(const T& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership
template <typename T> T* new_move(T&& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership

template <typename T> T* addr_if_possible(T& arg) {return &arg;}
template <typename T> T* addr_if_possible(T&& arg) {return nullptr;}

#define HEXED_QUAL_PTR_ACCESS(CONST, GET) \
  CONST T* get() CONST {GET(CONST)} \
  CONST T& operator*() CONST {return *get();} \
  CONST T* operator->() CONST {return get();} \
  CONST T& value() CONST \
  { \
    CONST T* data = get(); \
    HEXED_ASSERT(data, "pointer object is null"); \
    return *data; \
  } \

#define HEXED_PTR_ACCESS(GET) \
  HEXED_QUAL_PTR_ACCESS(, GET) \
  HEXED_QUAL_PTR_ACCESS(const, GET) \

}
#endif
