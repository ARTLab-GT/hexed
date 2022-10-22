#ifndef HEXED_KERNEL_FACTORY_HPP_
#define HEXED_KERNEL_FACTORY_HPP_

#include <memory>
#include "config.hpp"
#include "Vector_view.hpp"
#include "Stopwatch_tree.hpp"

/*
 * The purpose of this file is to provide a way to specify the number of dimensions and
 * row size of the kernels at run-time, whereas in the kernel implementations these are
 * compile-time parameters. The idea is that if we require that the parameters be within
 * a finite range, we can simply instantiate the kernel templates for every possible
 * combination of parameters. This is obviously not efficient in terms of compilation
 * work load or executable size, but it provides optimal execution speed and run-time
 * flexibility, which for our purpose is more important.
 */

namespace hexed
{

/*
 * Base class for a "kernel", a word which here means a callable object which accepts dimensionality
 * and row size as template arguments to improve performance.
 */
template <typename T, typename U = void>
class Kernel
{
  public:
  virtual ~Kernel() = default;
  virtual U operator()(Sequence<T>&) = 0;
  virtual U operator()(Sequence<T>& sequence, Stopwatch_tree& tree)
  {
    tree.work_units_completed += sequence.size();
    Stopwatch::Operator oper (tree.stopwatch);
    return (*this)(sequence);
  }
  virtual U operator()(Sequence<T>& sequence, Stopwatch_tree& category_tree, std::string name)
  {
    auto& specific_tree {category_tree.children.at(name)};
    specific_tree.work_units_completed += sequence.size();
    Stopwatch::Operator cat_oper (category_tree.stopwatch);
    Stopwatch::Operator spc_oper (specific_tree.stopwatch);
    return (*this)(sequence);
  }
  typedef std::unique_ptr<Kernel<T, U>> ptr_t; // default definition of `ptr_t` below
};

namespace kernel_lookup
{

// declares a pointer type which can point to any template instantiation.
// must be defined by each kernel class.
template <template<int, int> typename kernel>
using ptr_t = typename kernel<1, 2>::ptr_t;

/*
 * A recursive template to create the desired kernel instantiations. You can think
 * of it as a compile-time linked list, where each template contains another template
 * with one less row size. When the row size is 1, it wraps around to `config::max_row_size`
 * but with one less dimension via specialization below. The second specialization below
 * terminates the recursion when both the row size and the dimension are 1.
 */
template <template<int, int> typename kernel, int max_n_dim, int max_row_size>
class Kernel_lookup
{
  Kernel_lookup<kernel, max_n_dim, max_row_size - 1> decremented; // the next "link" in the recursive hirarchy
  public:
  // traverses the recursive hierarchy to find the desired template instance
  template<typename... constructor_args>
  ptr_t<kernel> get(int n_dim, int row_size, constructor_args&&... args)
  {
    if ((n_dim == max_n_dim) && (row_size == max_row_size)) return ptr_t<kernel>{new kernel<max_n_dim, max_row_size>(args...)}; // if I have the template you're looking for, return it
    else return decremented.get(n_dim, row_size, args...); // if I don't have it, delegate it to the next link (which is to say, "move along")
  }
};

// specialization to make the recursion wrap around to `config::max_row_size` of one less dimension
template <template<int, int> typename kernel, int max_n_dim>
class Kernel_lookup<kernel, max_n_dim, 1>
{
  Kernel_lookup<kernel, max_n_dim - 1, config::max_row_size> decremented;
  public:
  template<typename... constructor_args>
  ptr_t<kernel> get(int n_dim, int row_size, constructor_args&&... args)
  {
    if (n_dim == max_n_dim) return ptr_t<kernel>{new kernel<max_n_dim, 1>(args...)};
    else return decremented.get(n_dim, row_size, args...);
  }
};

// recursive base case
template <template<int, int> typename kernel>
class Kernel_lookup<kernel, 1, 1>
{
  public:
  template<typename... constructor_args>
  ptr_t<kernel> get(int n_dim, int row_size, constructor_args&&... args)
  {
    return ptr_t<kernel>{new kernel<1, 1>(args...)};
  }
};

} // namespace kernel_lookup

/*
 * Gets a `std::unique_ptr` to a `base_t` of `kernel`. The underlying kernel will be
 * instantiated with `n_dim` and `row_size` provided as template arguments. Arguments
 * must satisfy `3 >= n_dim > 0` and `config::max_row_size >= row_size > 1`.
 */
template <template<int, int> typename kernel, typename... constructor_args>
kernel_lookup::ptr_t<kernel> kernel_factory(int n_dim, int row_size, constructor_args&&... args)
{
  if ((n_dim < 1) || (n_dim > 3) || (row_size < 2) || (row_size > config::max_row_size)) {
    throw std::runtime_error("demand for invalid kernel");
  }
  kernel_lookup::Kernel_lookup<kernel, 3, config::max_row_size> lookup;
  return lookup.get(n_dim, row_size, args...);
}

}
#endif
