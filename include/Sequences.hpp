#ifndef HEXED_SEQUENCES_HPP_
#define HEXED_SEQUENCES_HPP_

#include "Sequence.hpp"

namespace hexed
{

template<typename S, typename T>
S trivial_convert(T& t)
{
  return t;
}

template <typename ref_t, typename ptr_t>
ref_t ptr_convert(ptr_t& ptr)
{
  return *ptr;
}

/*! \brief Provides a view of an existing container (e.g., `std::vector`) as a `Sequence`.
 * \details `Vector_view` also allows the vector elements to be viewed as a different type than they
 * are stored as with a user-supplied conversion function
 * (again, useful if you don't to provide access to only some aspects of the underlying data).
 * \note if write acces is desired, `reference_t` should be a reference type.
 */
template<typename reference_t, typename storage_t = reference_t,
         reference_t (*convert)(storage_t&) = &trivial_convert<reference_t, storage_t>,
         template<typename> typename sequence_template = std::vector>
class Vector_view : public Sequence<reference_t>
{
  sequence_template<storage_t>& vec;

  public:
  Vector_view(sequence_template<storage_t>& viewed)
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

/*! \brief A `Sequence` formed by concatenating two `Sequence`s.
 * \details For better or for worse, no copies are made---it requres references to existing sequences.
 * More then 2 sequences can be concatenated by nesting `Concatenation`s.
 */
template <typename T>
class Concatenation : public Sequence<T>
{
  Sequence<T>* s0;
  Sequence<T>* s1;

  public:
  Concatenation(Sequence<T>& seq0, Sequence<T>& seq1)
  : s0{&seq0}, s1{&seq1}
  {}

  int size()
  {
    return s0->size() + s1->size();
  }

  T operator[](int index)
  {
    int index1 = index - s0->size();
    return (index1 < 0) ? (*s0)[index] : (*s1)[index1];
  }
};

//! \brief Provides a view of some elements of a `Sequence` identified by a sequence of indices.
template <typename T, typename seq_t = std::vector<int>>
class Index : public Sequence<T>
{
  Sequence<T>& s;
  seq_t i;
  public:
  Index(Sequence<T>& seq, seq_t inds) : s{seq}, i{inds} {}
  int size()
  {
    return i.size();
  }

  T operator[](int index)
  {
    return s[i[index]];
  }
};

}
#endif
