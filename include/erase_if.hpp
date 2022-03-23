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
}

}
#endif
