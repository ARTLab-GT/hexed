#ifndef HEXED_MUTUAL_PTR_HPP_
#define HEXED_MUTUAL_PTR_HPP_

#include "assert.hpp"

namespace hexed
{

/*! \brief Abstract base class for implementing mutually-connected pointers.
 * \details Specifies an interface that has a reference to a `T` (`mine()`),
 * as well as a mechanism to connect with and disconnect from a `Ptr_base<T, U>`.
 * Implements the mechanics of mutual connection and disconnection.
 * Derived classes should implement `mine()`, `connect_self(Ptr_base<U, T>*)`, and `disconnect_self(Ptr_base<U, T>*)`
 * and then implement their own interface that calls `connect(Ptr_base<U, T>*)` and `disconnect(Ptr_base<U, T>*)`.
 */
template <typename T, typename U>
class Ptr_base
{
  friend class Ptr_base<U, T>;

  protected:
  //! \brief Do whatever I need to do to connect myself to this `Ptr_base` without worrying about the other side.
  virtual void connect_self(Ptr_base<U, T>*) = 0;
  //! \brief Do whatever I need to do to disconnect myself from this `Ptr_base` without worrying about the other side.
  //! \details What happends if I'm not already connected to the `Ptr_base` is up to the derived class.
  virtual void disconnect_self(Ptr_base<U, T>*) = 0;

  //! \brief Mutually connect us by calling both our `connect_self`.
  void connect(Ptr_base<U, T>* other)
  {
    connect_self(other);
    if (other) other->connect_self(this);
  }

  //! \brief Mutually disconnect us by calling both our `disconnect_self`.
  void disconnect(Ptr_base<U, T>* other)
  {
    if (other) other->disconnect_self(this);
    disconnect_self(other);
  }

  public:
  virtual ~Ptr_base() = default;
  virtual T& mine() = 0; //!< \brief Obtain the object I am permanently associated with.
  virtual const T& mine() const = 0; //!< \overload
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
 * and can point to objects of type `U` by pairing with a `Mutual_ptr<U, T>`, or more generally, any `Ptr_base<U, T>`.
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
  //! \see `operator=(Mutual_ptr&& other)`
  Mutual_ptr(Mutual_ptr&& other) : _mine{nullptr}, _partner{nullptr} {*this = std::move(other);}

  /*! \brief Move semantics steal the other's pairing.
   * \details If `other` is unpaired, `this` will be unpaired.
   * `other` is left unpaired with its `mine` intact.
   */
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

  #define ACCESS(CONST) \
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

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

/*! \brief Like a `Mutual_ptr`, except it can accept multiple partners.
 * \details Has a sequence of `partners()`, each of which is a `Ptr_base<U, T>`.
 * Thus it can be connected to `Mutual_ptr`s or other `Multiple_ptr`s.
 * Disconnecting (using `Multiple_ptr::remove`, `Mutual_ptr::unpair()`, or destroying the partners)
 * will simply remove partners from the sequence.
 */
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
  //! \brief Constructs a `Multiple_ptr` and permanently associates it with `data`.
  //! \details `data` must not be null.
  Multiple_ptr(T* data) : _mine{data} {HEXED_ASSERT(_mine, "`Multiple_ptr` cannot be constructed from null data.");}
  ~Multiple_ptr() {for (auto p : _partners) this->disconnect(p);}
  //! \brief Cannot copy a `Multiple_ptr` because it is unclear what that should do to partners.
  Multiple_ptr(const Multiple_ptr&) = delete;
  Multiple_ptr& operator=(const Multiple_ptr&) = delete;
  //! \see `operator=(Multiple_ptr&&)`
  Multiple_ptr(Multiple_ptr&& other) : _mine{nullptr} {*this = std::move(other);}

  /*! \brief Moving a `Multiple_ptr` steals the others partners.
   * \details `mine()` will point to `other.mine()`.
   * Any partners `this` had before the move assignment will be removed.
   * `other` will be left with its `mine()` intact and no partners.
   */
  Multiple_ptr& operator=(Multiple_ptr&& other)
  {
    _mine = other._mine;
    for (auto p : _partners) this->disconnect(p);
    for (auto p : other._partners) {
      other.disconnect(p);
      this->connect(p);
    }
    return *this;
  }

  //! \brief Adds `other` to the list of partners and reciprocally connects it with `this`.
  void add(Ptr_base<U, T>& other) {this->connect(&other);}
  //! \brief Removes `other` from the list of partners and reciprocally disconnects it from `this`.
  void remove(Ptr_base<U, T>& other) {this->disconnect(&other);}
  //! \brief Returns `true` iff `this` has at least 1 partner.
  operator bool() const {return !_partners.empty();}

  #define ACCESS(CONST) \
    CONST T& mine() CONST {return *_mine;} \
    /*! \brief Obtains the sequence of partners, in no particular order. */ \
    std::vector<CONST Ptr_base<U, T>*> partners() CONST { \
      return _partners; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
