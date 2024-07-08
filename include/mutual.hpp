#ifndef HEXED_MUTUAL_HPP_
#define HEXED_MUTUAL_HPP_

#include <vector>
#include "assert.hpp"
#include "Sequence.hpp"

/*! \brief a namespace for mutually-connected objects
 * \details When implementing adaptive meshing, there are many situations where:
 * - Two or more objects need to know about each other.
 * - One of the objects may be moved or destroyed,
 *   in which case the other objects need to be notified so that they are not left with dangling pointers.
 *
 * `std::shared_ptr` does not solve this problem.
 * I am not trying to _preserve_ an object until there are no remaining references to it.
 * When it is time for an object to be deleted, I simply want everyone currently referencing it to _know_.
 * This namespace provides the fundamental base classes to meet this need.
 */
namespace hexed::mutual {

/*! \brief abstract base class for all mutually-connected objects
 * \details Implements the basic mechanics of mutual connection and disconnection.
 * Derived classes may provide data of currently unknown type,
 * so `Base` takes a template parameter `T` and provides a virtual function `_mine()`
 * to access data of type `T`.
 * It may be desireable for a `Base` to be mutually connected with another `Base` (its "partner")
 * that provides data of a different type, so `Base` takes a second parameter `U` which is its partner's `T`.
 * Derived classes must override `_set()` and `_unset()` to implement the process of connecting themself to a partner.
 * Derived classes should implement an interface that calls `_connect` and `_disconnect` as desired, including:
 * - A destructor that `_disconnect`s the object from its partners.
 * - Move semantics that cause the moved-to object to steal the partners of the moved-from object.
 *
 * Derived classes may override `_mine()` to provide their partners access to some data.
 */
template <typename T, typename U>
class Base {
  friend class Base<U, T>;
protected:
  //! \brief takes the necessary steps to connect `this` to `that`, without worrying about anything on `that`'s end
  virtual void _set(Base<U, T>& that) = 0;
  //! \brief takes the necessary steps to disconnect `this` from `that`, without worrying about anything on `that`'s end
  virtual void _unset(Base<U, T>& that) = 0;
  //! \brief may be overridden by derrived classes to provide partners to data of some arbitrary type `T`
  virtual T* _mine() {return nullptr;}
  virtual const T* _mine() const {return nullptr;} //!< \overload

  //! \brief mutually connects `this` and `that` by calling both of their `_set()` member functions
  void _connect(Base<U, T>& that) {
    _set(that);
    that._set(*this);
  }

  //! \brief mutually disconnects `this` and `that` by calling both of their `_unset()` member functions
  void _disconnect(Base<U, T>& that) {
    that._unset(*this);
    _unset(that);
  }

  //! \brief Accesses the `_mine()` of `that`
  static U* _yours(Base<U, T>& that) {return that._mine();}
  static const U* _yours(const Base<U, T>& that) {return that._mine();} //!< \overload
};

//! \brief an object which is mutually connected ("paired") with only one other object
template <typename T, typename U>
class Single : public Base<T, U> {
public:
  //! \brief Constructs a `Single` with no partner (it is "unpaired")
  Single() : _partner{nullptr} {}
  Single(const Single&) = delete;
  //! \brief steals `that`'s partner, if it has one, leaving `that` unpaired
  Single(Single&& that) : _partner{nullptr} {*this = std::move(that);}
  ~Single() {unpair();}
  Single& operator=(const Single&) = delete;

  //! \brief disconnects `this` from its partner, if it has one, and steals `that`'s, if it has one
  //! \details leaves `that` unpaired
  Single& operator=(Single&& that) {
    unpair();
    if (that._partner) {
      pair(*that._partner);
      that.unpair();
    }
    return *this;
  }

  //! \brief mutually connects `this` with `that`
  void pair(Base<U, T>& that) {this->_connect(that);}
  //! \brief mutually disconnects `this` from `that`
  void unpair() {if (_partner) this->_disconnect(*_partner);}
  //! \brief returns `true` iff this is currently paired with a partner
  bool paired() const {return _partner;}
  //! \brief provides access to `this`'s partner
  Base<U, T>* partner() {return _partner;}
  const Base<U, T>* partner() const {return _partner;} //!< \overload

private:
  void _set(Base<U, T>& other) override {
    unpair();
    _partner = &other;
  }
  void _unset(Base<U, T>& other) override {
    HEXED_ASSERT(_partner == &other, "not connected to `other`");
    _partner = nullptr;
  }
  Base<U, T>* _partner;
};

/*! \brief an object that can be mutually connected with any number of partners
 * \details This class does not provide any interface for connecting to and disconnecting from partners.
 * It is intended that you either connect `Single`'s to it,
 * or create a derrived class which provides an interface to connect and disconnect
 * (by calling `Base::_connect` and `Base::_disconnect`).
 * It does, however, take care of destruction and move semantics,
 * so derrived classes shouldn't have to worry about that.
 */
template <typename T, typename U>
class Multiple : public Base<T, U>
{
public:
  Multiple() = default;
  Multiple(const Multiple&) = delete;
  //! \brief steals all of `that`'s partners, leaving `that` unconnected
  Multiple(Multiple&& that) {*this = std::move(that);}
  ~Multiple() {for (int i = _partners.size() - 1; i >= 0; --i) this->_disconnect(*_partners[i]);}
  Multiple& operator=(const Multiple&) = delete;

  //! \brief steals all of `that`'s partners, leaving `that` unconnected
  //! \details `this` is first disconnected from all its own partners
  Multiple& operator=(Multiple&& that) {
    for (int i = _partners.size() - 1; i >= 0; --i) this->_disconnect(*_partners[i]);
    for (int i = that._partners.size() - 1; i >= 0; --i) {
      Base<U, T>* p = that._partners[i];
      that._disconnect(*p);
      this->_connect(*p);
    }
    return *this;
  }

  #define ACCESS(CONST) \
    /*! \brief provides access to the list of partners */ \
    next::Sequence<CONST Base<U, T>&> partners() CONST { \
      return { \
        [this](std::size_t index)->CONST Base<U, T>& {return *_partners[index];}, \
        [this](){return _partners.size();}, \
      }; \
    }
  ACCESS()
  ACCESS(const)
  #undef ACCESS
  std::vector<Base<U, T>*> _partners;

private:
  void _set(Base<U, T>& that) override {
    if (!std::any_of(_partners.begin(), _partners.end(), [&that](Base<U, T>* p){return p == &that;})) {
      _partners.push_back(&that);
    }
  }
  void _unset(Base<U, T>& that) override {std::erase(_partners, &that);}
};

}
#endif
