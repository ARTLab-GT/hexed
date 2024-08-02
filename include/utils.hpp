#ifndef HEXED_UTILS_HPP_
#define HEXED_UTILS_HPP_

#include <iostream>
#include "assert.hpp"

namespace hexed {

template <typename T> T gotcha(T) {throw std::runtime_error("hexed::gotcha");} //!< \brief don't ask

//! \brief prints its argument and then returns it
//! \details you can wrap this around an expression to print it without computing it again
template <typename T>
T& printed(T& t) {
  std::cout << t << std::endl;
  return t;
}

template <typename T> T* new_copy(const T& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership
template <typename T> T* new_move(T&& t) {return new T(t);} //!< useful for cppyy which doesn't like to relinquish ownership

template <typename T> std::add_pointer<T>::type addr_if_possible(T& arg) {return &arg;}
template <typename T> std::add_pointer<T>::type addr_if_possible(T&& arg) {
  HEXED_THROW(format_str(200, "attempt to take address of non-addressable type %s", typeid(T).name()));
}

//! \brief gets the extension of a file name and converts it to lowercase
std::string file_extension(std::string file_name);

}
#endif
