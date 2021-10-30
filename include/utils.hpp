#ifndef CARTDG_UTILS_HPP_
#define CARTDG_UTILS_HPP_

#include <iostream>
#include <string>

namespace cartdg
{

template <typename T>
T print(T value, std::string prefix="", std::string suffix="\n")
{
  std::cout << prefix << value << suffix;
  return value;
}

}
#endif
