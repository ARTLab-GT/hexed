#ifndef HEXED_MUTUAL_PTR_HPP_
#define HEXED_MUTUAL_PTR_HPP_

#include "assert.hpp"

namespace hexed
{

template <typename T, typename U>
class Ptr_base
{
  friend class Ptr_base<U, T>;
  protected:
  virtual void connect_self(Ptr_base<U, T>*) = 0;
  virtual void disconnect_self(Ptr_base<U, T>*) = 0;
  void connect(Ptr_base<U, T>* other)
  {
    connect_self(other);
    if (other) other->connect_self(this);
  }
  void disconnect(Ptr_base<U, T>* other)
  {
    if (other) other->disconnect_self(this);
    disconnect_self(other);
  }
  public:
  virtual ~Ptr_base() = default;
  virtual T& mine() = 0;
  virtual const T& mine() const = 0;
};

/*! \brief For creating pairs of mutually connected pointers in a robust way.
 * \details There are several instances in Hexed where two objects need to reference each other.
 * If either one of them is deleted, the other needs to be notified to avoid creating a dangling pointer.
 * This class implements that idea by creating a type of pointer
 * that can be paired with another pointer to create a mutual reference,
 * but breaking the connection from either end will mutually disconnect both pointers.
 * A `Mutual_ptr` does not own any data.
 * It only manages pairing and unpairing.
 * A `Mutual_ptr<T, U>` is permanently associated with an object of type `T`
 * and can point to objects of type `U` by pairing with a `Mutual_ptr<U, T>`.
 */
template <typename T, typename U>
class Mutual_ptr : public Ptr_base<T, U>
{
  T* _mine;
  Ptr_base<U, T>* _partner;
  void connect_self(Ptr_base<U, T>* other) override
  {
    unpair();
    _partner = other;
  }
  void disconnect_self(Ptr_base<U, T>* other) override
  {
    HEXED_ASSERT(_partner == other, "can only disconnect a `Mutual_ptr` from its current partner");
    _partner = nullptr;
  }

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
  void pair(Ptr_base<U, T>& other) {this->connect(&other);}

  //! \brief if `this` is in a pair, mutually unpairs `this` and its partner
  //! \details safe to call regardless of whether `this` is currently paired
  void unpair() {this->disconnect(_partner);}

  operator bool() const {return _partner;} //!< \brief returns `true` iff currently paired
  #define ACCESS \
    /*! \brief the object that `this` is permanently associated with */ \
    CONST T& mine() CONST {return *_mine;} \
    /*! \brief if paired, the `Mutual_ptr` this is currently paired with; else `nullptr` */ \
    CONST Ptr_base<U, T>* partner() CONST {return _partner;} \
    /*! \brief if paired, returns `partner()`'s `mine()`; else `nullptr` */ \
    CONST U* get() CONST {return _partner ? &_partner->mine() : nullptr;} \
    /*! \brief dereferencing obtains `partner()`s `mine()` (undefined if unpaired) */ \
    CONST U& operator*() CONST {return _partner->mine();} \
    /*! \brief dereferencing members refers to `partner()`s `mine()` (undefined if unpaired) */ \
    CONST U* operator->() CONST {return &_partner->mine();} \
    /*! \brief obtains reference to `partner()`'s `mine()` if paired and throws exception if not paired */ \
    CONST U& value()  CONST \
    { \
      HEXED_ASSERT(partner(), "attempt to get the value of an unpaired `Mutual_ptr`"); \
      return _partner->mine(); \
    } \

  #define CONST
  ACCESS
  #undef CONST
  #define CONST const
  ACCESS
  #undef CONST
  #undef ACCESS
};

template <typename T, typename U>
class Multiple_ptr : public Ptr_base<T, U>
{
  T* _mine;
  std::vector<Ptr_base<U, T>*> _partners;
  void connect_self(Ptr_base<U, T>* other) override
  {
    _partners.push_back(other);
  }
  void disconnect_self(Ptr_base<U, T>* other) override
  {
    std::erase(_partners, other);
  }

  public:
  Multiple_ptr(T* data) : _mine{data} {HEXED_ASSERT(_mine, "`Multiple_ptr` cannot be constructed from null data.");}
  ~Multiple_ptr() {for (auto p : _partners) this->disconnect(p);}
  void add(Ptr_base<U, T>& other) {this->connect(&other);}
  void remove(Ptr_base<U, T>& other) {this->disconnect(&other);}

  #define ACCESS \
    CONST T& mine() CONST {return *_mine;} \
    std::vector<CONST Ptr_base<U, T>*> partners() CONST { \
      return _partners; \
    } \

  #define CONST
  ACCESS
  #undef CONST
  #define CONST const
  ACCESS
  #undef CONST
  #undef ACCESS
};

}
#endif
