#ifndef HEXED_MUTUAL_PTR_HPP_
#define HEXED_MUTUAL_PTR_HPP_

#include "assert.hpp"

namespace hexed
{

template <typename T, typename U>
class Mutual_ptr
{
  T* _mine;
  Mutual_ptr<U, T>* _partner;
  friend Mutual_ptr<U, T>;

  public:
  Mutual_ptr(T* data) : _mine{data}, _partner{nullptr}
  {
    HEXED_ASSERT(_mine, "`Mutual_ptr` cannot be constructed from null data");
  }
  Mutual_ptr(const Mutual_ptr&) = delete;
  Mutual_ptr(Mutual_ptr&& other) {*this = other;}
  ~Mutual_ptr() {unpair();}
  Mutual_ptr& operator=(const Mutual_ptr&) = delete;
  Mutual_ptr& operator=(Mutual_ptr&& other)
  {
    _mine = other._mine;
    _partner = nullptr;
    if (other._partner) {
      auto p = other._partner;
      p->unpair();
      pair(*p);
    }
    return *this;
  }

  void pair(Mutual_ptr<U, T>& other)
  {
    unpair();
    other.unpair();
    _partner = &other;
    other._partner = this;
  }
  void unpair()
  {
    if (_partner) {
      _partner->_partner = nullptr;
      _partner = nullptr;
    }
  }

  operator bool() {return _partner;}
  T* mine() {return _mine;}
  Mutual_ptr<U, T>* partner() {return _partner;}
  U* get() {return _partner ? _partner->_mine : nullptr;}
  operator U*() {return get();}
  U& operator*() {return *_partner->_mine;}
  U* operator->() {return _partner->_mine;}
  U& value()
  {
    HEXED_ASSERT(partner(), "attempt to get the value of an unpaired `Mutual_ptr`");
    return *this;
  }
};

}
#endif
