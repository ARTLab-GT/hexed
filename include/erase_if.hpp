#ifndef HEXED_ERASE_IF_HPP_
#define HEXED_ERASE_IF_HPP_

#include <vector>

namespace hexed
{

// reimplementation of `std::erase_if` which we can't use
// because we don't have c++20 on the lab machines
template <typename T, typename callable>
void erase_if(std::vector<T>& vec, callable condition)
{
  auto new_end = std::remove_if(vec.begin(), vec.end(), condition);
  vec.erase(new_end, vec.end());
}

}
#endif
