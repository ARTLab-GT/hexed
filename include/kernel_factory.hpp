#ifndef CARTDG_KERNEL_FACTORY_HPP_
#define CARTDG_KERNEL_FACTORY_HPP_

#include <memory>
#include <cartdgConfig.hpp>

namespace cartdg
{

template <template<int, int, int> typename kernel>
class Kernel_traits
{
  public:
  class base_t;
};

namespace kernel_lookup
{

template <template<int, int, int> typename kernel>
using ptr_t = std::unique_ptr<typename Kernel_traits<kernel>::base_t>;

template <template<int, int, int> typename kernel, int n_dim, int row_size>
ptr_t<kernel> pointer()
{
  return ptr_t<kernel>{new kernel<n_dim+2, n_dim, row_size>};
}

template <template<int, int, int> typename kernel, int max_n_dim, int max_row_size>
class Kernel_lookup
{
  Kernel_lookup<kernel, max_n_dim, max_row_size - 1> decremented;
  public:
  ptr_t<kernel> get(int n_dim, int row_size)
  {
    if ((n_dim == max_n_dim) && (row_size == max_row_size)) return pointer<kernel, max_n_dim, max_row_size>();
    else return decremented.get(n_dim, row_size);
  }
};

template <template<int, int, int> typename kernel, int max_n_dim>
class Kernel_lookup<kernel, max_n_dim, 1>
{
  Kernel_lookup<kernel, max_n_dim - 1, config::max_row_size> decremented;
  public:
  ptr_t<kernel> get(int n_dim, int row_size)
  {
    if (n_dim == max_n_dim) return pointer<kernel, max_n_dim, 1>();
    else return decremented.get(n_dim, row_size);
  }
};

template <template<int, int, int> typename kernel>
class Kernel_lookup<kernel, 1, 1>
{
  public:
  ptr_t<kernel> get(int n_dim, int row_size)
  {
    return pointer<kernel, 1, 1>();
  }
};

} // namespace kernel_lookup

template <template<int, int, int> typename kernel>
kernel_lookup::ptr_t<kernel> kernel_factory(int n_dim, int row_size)
{
  kernel_lookup::Kernel_lookup<kernel, 3, config::max_row_size> lookup;
  return lookup.get(n_dim, row_size);
}

}
#endif
