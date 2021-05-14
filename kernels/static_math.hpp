#ifndef CARTDG_STATIC_MATH_HPP_
#define CARTDG_STATIC_MATH_HPP_

namespace cartdg
{

namespace static_math
{

template<typename number_t>
constexpr number_t pow(number_t base, int exponent)
{
  number_t result = 1;
  for (int i = 0; i < exponent; ++i) result *= base;
  for (int i = 0; i > exponent; --i) result /= base;
  return result;
}

}
}
#endif
