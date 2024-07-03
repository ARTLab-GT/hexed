#ifndef HEXED_MUTUAL_HPP_
#define HEXED_MUTUAL_HPP_

#include <vector>
#include "assert.hpp"

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

  void _connect(Base<U, T>& other)
  {
    _set(other);
    other._set(*this);
  }

  void _disconnect(Base<U, T>& other)
  {
    other._unset(*this);
    _unset(other);
  }
};

template <typename T, typename U>
class Single : public Base<T, U>
{
  Base<U, T>* _partner;

  void _set(Base<U, T>& other) override
  {
    unpair();
    _partner = &other;
  }

  void _unset(Base<U, T>& other) override
  {
    HEXED_ASSERT(_partner == &other, "not connected to `other`");
    _partner = nullptr;
  }

  protected:
  Base<U, T>* _get() {return _partner;}
  const Base<U, T>* _get() const {return _partner;}

  public:
  Single() : _partner{nullptr} {}
  Single(const Single&) = delete;
  Single(Single&& other) : _partner{nullptr} {*this = std::move(other);}
  ~Single() {unpair();}
  Single& operator=(const Single&) = delete;

  Single& operator=(Single&& other)
  {
    unpair();
    if (other._partner) {
      pair(*other._partner);
      other.unpair();
    }
    return *this;
  }

  void pair(Base<U, T>& other) {this->_connect(other);}
  void unpair() {if (_partner) this->_disconnect(*_partner);}
};

template <typename T, typename U>
class Multiple : public Base<T, U>
{
  std::vector<Base<U, T>*> _partners;

  void _set(Base<U, T>& other) override {_partners.push_back(&other);}
  void _unset(Base<U, T>& other) override {std::erase(_partners, &other);}

  protected:
  std::vector<Base<U, T>*> _get() {return _partners;}
  const std::vector<Base<U, T>*> _get() const {return _partners;}

  public:
  Multiple() = default;
  Multiple(const Multiple&) = delete;
  Multiple(Multiple&& other) {*this = std::move(other);}
  ~Multiple() {for (auto p : _partners) this->_disconnect(*p);}
  Multiple& operator=(const Multiple&) = delete;

  Multiple& operator=(Multiple&& other)
  {
    for (auto p : _partners) this->_disconnect(*p);
    for (auto p : other._partners) {
      other._disconnect(*p);
      this->_connect(*p);
    }
    return *this;
  }
};

}
#endif
