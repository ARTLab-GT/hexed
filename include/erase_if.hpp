#ifndef CARTDG_ERASE_IF_HPP_
#define CARTDG_ERASE_IF_HPP_

#include <vector>

namespace cartdg
{

/*
 * Erase every element `x` in a vector where `condition(x)` evaluates to `true`.
 * Complexity is O(`vec.size()`/(number of threads) + (number of threads))
 * regardless of number of elements to be deleted.
 * Thus this is much more efficient than repeatedly calling `vec.erase` if
 * multiple elements which are not in a contiguous range are to be erased.
 */
template <typename T>
void erase_if(std::vector<T>& vec, bool (*condition)(const T&))
{
  // FIXME: make this parallel
  std::vector<T> erased;
  for (unsigned i_elem = 0; i_elem < vec.size(); ++i_elem) {
    if (!condition(vec[i_elem])) erased.push_back(std::move(vec[i_elem]));
  }
  vec = std::move(erased);
}

}
#endif
