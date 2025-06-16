#ifndef HEXED_FORMAT_STR_HPP_
#define HEXED_FORMAT_STR_HPP_

namespace hexed {

/*! \brief Standard string formatting.
 * \details Basically a knockoff of `std::format` in C++20 (which at the time of writing we can't use on the lab machines).
 * Invokes `snprintf`, but works in terms of `std::string`s and handles buffer creation for you.
 * If length of resulting string is greater than `max_chars`, it will iteratively allocate more characters.
 */
template <typename... format_args>
std::string format_str(int max_chars, std::string fstring, format_args... args) {
  std::vector<char> buffer(max_chars);
  int overflow;
  do {
    overflow = snprintf(buffer.data(), max_chars, fstring.c_str(), args...);
    if (overflow < 0) throw std::runtime_error("encoding error in `hexed::format_str`");
    max_chars *= 2;
    buffer.resize(max_chars);
  } while (overflow >= max_chars);
  return std::string(buffer.data());
}

//! \brief sets `max_chars` to 200.
template <typename... format_args>
std::string format_str(std::string fstring, format_args... args) {
  return format_str(200, fstring, args...);
}

}
#endif
