#ifndef HEXED_MUTUAL_HPP_
#define HEXED_MUTUAL_HPP_

#include <vector>
#include "assert.hpp"
#include "Sequence.hpp"

namespace hexed::mutual
{

template <typename T, typename U>
class Base
{
  friend class Base<U, T>;

  protected:
  virtual void _set(Base<U, T>&) = 0;
  virtual void _unset(Base<U, T>&) = 0;
  virtual T* _mine() {return nullptr;}
  virtual const T* _mine() const {return nullptr;}

  void _connect(Base<U, T>& other) {
    _set(other);
    other._set(*this);
  }

  void _disconnect(Base<U, T>& other) {
    other._unset(*this);
    _unset(other);
  }

  static U* _yours(Base<U, T>& other) {return other._mine();}
  static const U* _yours(const Base<U, T>& other) {return other._mine();}
};

template <typename T, typename U>
class Single : public Base<T, U>
{
  Base<U, T>* _partner;

  void _set(Base<U, T>& other) override {
    unpair();
    _partner = &other;
  }

  void _unset(Base<U, T>& other) override {
    HEXED_ASSERT(_partner == &other, "not connected to `other`");
    _partner = nullptr;
  }

  public:
  Single() : _partner{nullptr} {}
  Single(const Single&) = delete;
  Single(Single&& other) : _partner{nullptr} {*this = std::move(other);}
  ~Single() {unpair();}
  Single& operator=(const Single&) = delete;

  Single& operator=(Single&& other) {
    unpair();
    if (other._partner) {
      pair(*other._partner);
      other.unpair();
    }
    return *this;
  }

  void pair(Base<U, T>& other) {this->_connect(other);}
  void unpair() {if (_partner) this->_disconnect(*_partner);}
  bool paired() const {return _partner;}
  Base<U, T>* partner() {return _partner;}
  const Base<U, T>* partner() const {return _partner;}
};

template <typename T, typename U>
class Multiple : public Base<T, U>
{
  std::vector<Base<U, T>*> _partners;

  void _set(Base<U, T>& other) override {
    if (!std::any_of(_partners.begin(), _partners.end(), [&other](Base<U, T>* p){return p == &other;})) {
      _partners.push_back(&other);
    }
  }

  void _unset(Base<U, T>& other) override {std::erase(_partners, &other);}

  public:
  Multiple() = default;
  Multiple(const Multiple&) = delete;
  Multiple(Multiple&& other) {*this = std::move(other);}
  ~Multiple() {for (int i = _partners.size() - 1; i >= 0; --i) this->_disconnect(*_partners[i]);}
  Multiple& operator=(const Multiple&) = delete;

  Multiple& operator=(Multiple&& other) {
    for (int i = _partners.size() - 1; i >= 0; --i) this->_disconnect(*_partners[i]);
    for (int i = other._partners.size() - 1; i >= 0; --i) {
      Base<U, T>* p = other._partners[i];
      other._disconnect(*p);
      this->_connect(*p);
    }
    return *this;
  }

  #define ACCESS(CONST) \
    next::Sequence<CONST Base<U, T>&> partners() CONST { \
      return { \
        [this](std::size_t index)->CONST Base<U, T>& {return *_partners[index];}, \
        [this](){return _partners.size();}, \
      }; \
    }
  ACCESS()
  ACCESS(const)
  #undef ACCESS
};

}
#endif
