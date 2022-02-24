#ifndef CARTDG_VECTOR_VIEW_HPP_
#define CARTDG_VECTOR_VIEW_HPP_

namespace cartdg
{

template<typename S, typename T>
S trivial_convert(T& t)
{
  return t;
}

/*
 * Provides convenient but limited access to a `std::vector`. Supports read/write access
 * to the elements and evaluation of the size, but does not allow insertion or removal of
 * elements. This is useful for classes providing limited public access to private variables.
 * `Vector_view` also allows the vector elements to be viewed as a different type than they
 * are stored as with a user-supplied conversion function (again, useful if you don't to provide
 * access to only some aspects of the underlying data). Note: if write acces is desired,
 * `reference_t` should be a reference type.
 */
template<typename reference_t, typename storage_t = reference_t, reference_t (*convert)(storage_t&) = &trivial_convert<reference_t, storage_t>>
class Vector_view
{
  std::vector<storage_t>& vec;

  public:
  Vector_view(std::vector<storage_t>& viewed)
  : vec{viewed}
  {}

  int size()
  {
    return vec.size();
  }

  reference_t operator[](int index)
  {
    return convert(vec[index]);
  }
};

}
#endif
