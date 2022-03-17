#ifndef CARTDG_KERNEL_FACTORY_HPP_
#define CARTDG_KERNEL_FACTORY_HPP_

#include <memory>
#include <array>

namespace cartdg
{

class Base
{
  public:
  virtual int foo() = 0;
};

const int n = 5;

template <template<int> typename kernel, int i>
class Kernel_lookup
{
  Kernel_lookup<kernel, i - 1> decremented;
  public:
  std::unique_ptr<Base> get(int j)
  {
    if (i == j) return std::unique_ptr<Base>{new kernel<i>};
    else return decremented.get(j);
  }
};

template <template<int> typename kernel>
class Kernel_lookup<kernel, 0>
{
  public:
  std::unique_ptr<Base> get(int j)
  {
    return std::unique_ptr<Base>{new kernel<0>};
  }
};

template <template<int> typename kernel>
std::unique_ptr<Base> kernel_factory(int i)
{
  Kernel_lookup<kernel, n> lookup;
  return lookup.get(i);
}

template <int i>
class Derived : public Base
{
  public:
  virtual int foo() {return i;}
};

}
#endif
