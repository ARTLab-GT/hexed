#ifndef HEXED_ERASE_IF_HPP_
#define HEXED_ERASE_IF_HPP_

#include <vector>

namespace hexed
{

/*! \brief Erase every element `x` in a vector where `condition(x)` evaluates to `true`.
 * \details Complexity is O(`vec.size()`) regardless of number of elements to be deleted.
 * Thus this is much more efficient than repeatedly calling `vec.erase` if
 * multiple elements which are not in a contiguous range are to be erased.
 */
template <typename T, typename C>
void erase_if(std::vector<T>& vec, C condition)
{
  std::vector<T> erased;
  for (unsigned i_elem = 0; i_elem < vec.size(); ++i_elem) {
    if (!condition(vec[i_elem])) erased.push_back(std::move(vec[i_elem]));
  }
  vec = std::move(erased);
}

}
#endif
