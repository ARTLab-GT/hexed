#ifndef HEXED_MUTUAL_PTR_HPP_
#define HEXED_MUTUAL_PTR_HPP_

#include "assert.hpp"

namespace hexed
{

/*! \brief For creating pairs of mutually connected pointers in a robust way.
 * \details There are several instances in Hexed where two objects need to reference each other.
 * If either one of them is deleted, the other needs to be notified to avoid creating a dangling pointer.
 * This class implements that idea by creating a type of pointer that can be paired with another pointer to create a mutual reference,
 * but breaking the connection from either end will mutually disconnect both pointers.
 * A `Mutual_ptr` does not own any data.
 * It only manages pairing and unpairing.
 * A `Mutual_ptr<T, U>` is permanently associated with an object of type `T`
 * and can point to objects of type `U` by pairing with a `Mutual_ptr<U, T>`.
 */
template <typename T, typename U>
class Mutual_ptr
{
  T* _mine;
  Mutual_ptr<U, T>* _partner;
  friend Mutual_ptr<U, T>;

  public:
  //! \brief Constructs a `Mutual_ptr` permanently associated with `data`.
  //! \details `data` cannot be null. Henceforth, the value of `*data` can be changed but its address cannot.
  Mutual_ptr(T* data) : _mine{data}, _partner{nullptr}
  {
    HEXED_ASSERT(_mine, "`Mutual_ptr` cannot be constructed from null data");
  }
  ~Mutual_ptr() {unpair();} //!< \brief \ref unpair() "unpairs"

  //! \brief Copying is not supported, since what that ought to do to pairs is unclear.
  Mutual_ptr(const Mutual_ptr&) = delete;
  Mutual_ptr& operator=(const Mutual_ptr&) = delete;
  //! \brief Move semantics steal the other's pairing, if applicable
  Mutual_ptr(Mutual_ptr&& other) : _mine{nullptr}, _partner{nullptr}
  {
    *this = std::move(other);
  }
  Mutual_ptr& operator=(Mutual_ptr&& other)
  {
    unpair();
    _mine = other._mine;
    if (other._partner) pair(*other._partner);
    return *this;
  }

  //! \brief pair with another `Mutual_ptr`, breaking any existing pairs involving `this` or `other`
  //! \details safe to call regardless of whether `this` and `other` were previously paired
  void pair(Mutual_ptr<U, T>& other)
  {
    unpair();
    other.unpair();
    _partner = &other;
    other._partner = this;
  }
  //! \brief if `this` is in a pair, mutually unpairs `this` and its partner
  //! \details safe to call regardless of whether `this` is currently paired
  void unpair()
  {
    if (_partner) {
      _partner->_partner = nullptr;
      _partner = nullptr;
    }
  }

  operator bool() {return _partner;} //!< \brief returns `true` iff currently paired
  T& mine() {return *_mine;} //!< \brief the object that `this` is permanently associated with
  Mutual_ptr<U, T>* partner() {return _partner;} //!< \brief if paired, the `Mutual_ptr` this is currently paired with; else `nullptr`
  U* get() {return _partner ? _partner->_mine : nullptr;} //!< \brief if paired, returns `partner()`'s `mine()`; else `nullptr`
  operator U*() {return get();} //!< \brief casting `this` to `U*` returns the same as `get()`
  U& operator*() {return *_partner->_mine;} //!< \brief dereferencing obtains `partner()`s `mine()` (undefined if unpaired)
  U* operator->() {return _partner->_mine;} //!< \brief dereferencing members refers to `partner()`s `mine()` (undefined if unpaired)
  //! \brief obtains reference to `partner()`'s `mine()` if paired and throws exception if not paired
  U& value()
  {
    HEXED_ASSERT(partner(), "attempt to get the value of an unpaired `Mutual_ptr`");
    return *_partner->_mine;
  }
};

}
#endif
