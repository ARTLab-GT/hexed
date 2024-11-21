#ifndef HEXED_RECIPROCAL_HPP_
#define HEXED_RECIPROCAL_HPP_

#include "Mortal.hpp"

namespace hexed {

/*! \brief Implements pairs of reciprocally-connected pointers.
 * \details This class addresses the case where object `a` needs to have a pointer to object `b`
 * iff object `b` has a pointer to `a`.
 * This might sound like exactly what `mutual::Single` does,
 * but the problem is that a `Single<T, U>` can _only_ point to a `Single<U, T>`.
 * If you want an object to have a variety of reciprocal pointers to objects of multiple different types,
 * you need another layer of indirection:
 * `a` has a dedicated member which is a pointer object that derives from `Single` to point to `b`,
 * and `b` needs to have a pointer object that points to `a`,
 * and these pointers need to know about each other to enforce reciprocity.
 * Furthermore, both `a` and `b` must be derived from `Mortal`.
 * That is what this class provides.
 * If you want objects of class `A` and of class `B` to be able to point to each other,
 * give them members of type `Reciprocal_ptr<A, B>` and `Reciprocal_ptr<B, A>`, respectively,
 * and `pair()` the `Reciprocal_ptr`s.
 */
template <typename T, typename U>
class Reciprocal_ptr : public mutual::Single<T, U> {
public:
  //! \brief constructs a `Reciprocal_ptr` and sets its `mine` to `data` (which can be null)
  Reciprocal_ptr(T* data) : mine(data) {}
  operator bool() const {return get() != nullptr;} //!< \brief returns `true` iff `this` is not null

  #define ACCESS(CONST) \
    /*! \brief points to the `mine` of the `partner()` of `this` */ \
    /*! \details if not paired, returns `nullptr` */ \
    CONST U* get() CONST { \
      if (this->paired()) return this->_yours(*this->partner()); \
      return nullptr; \
    } \
    /*! \brief obtains a reference to the object `this` points to */ \
    /*!  \details undefined behavior if `this` is null */ \
    CONST U& operator*() CONST {return *get();} \
    /*! \brief accesses the members of the object `this` points to */ \
    /*! \details undefined behavior if `this` is null */ \
    CONST U* operator->() CONST {return get();} \
    /*! \brief obtains a reference to the object `this` points to */ \
    /*!  \details throws an exception if `this` is null */ \
    CONST U& value() CONST { \
      CONST U* data = get(); \
      HEXED_ASSERT(data, "`get()` is null"); \
      return *data; \
    } \

  ACCESS()
  ACCESS(const)
  #undef ACCESS

  /*! \brief The object that `this` is a proxy for in reciprocal connections.
   * \details Usually, `this` will be a member of the object `mine` points to,
   * but that is not technically necessary.
   * The user should feel free to modify this data member at will (including setting it to null).
   */
  Mortal_ptr<T> mine;

private:
  T* _mine() override {return mine.get();}
  const T* _mine() const override {return mine.get();}
};

/*! \brief Like `Reciprocal_ptr`, put it can be connected with multiple partners
 * \details `Reciprocal_list`s can be connected with other `Reciprocal_list`s or with `Reciprocal_ptr`s.
 * Of course, in both cases connections are reciprocated.
 */
template <typename T, typename U>
class Reciprocal_list : public mutual::Multiple<T, U> {
  public:
  //! \brief Constructs a `Reciprocal_list` and sets its `mine` to `data` (which can be null).
  Reciprocal_list(T* data) : mine(data) {}
  //! \brief Reciprocally connects `this` with `that`.
  //! \details In other words, adds `that` to `this`'s list of partners.
  void add(mutual::Base<U, T>& that) {this->_connect(that);}
  //! \brief disconnects all partners
  void clear() {for (auto& p : this->partners()) this->_disconnect(p);}

  //! \brief Reciprocally disconnects `this` from `that`.
  //! \details In other words, removes `that` from `this`'s list of partners.
  void remove(mutual::Base<U, T>& that) {
    if (this->partners().address().contains(&that)) this->_disconnect(that);
  }

  #define ACCESS(CONST) \
    /*! \brief Obtains all of `this`'s partners' `mine`s as a `next::Sequence` */ \
    next::Sequence<CONST U*> theirs() CONST { \
      auto p = this->partners(); \
      return { \
        [p](std::size_t index){return mutual::Base<T, U>::_yours(p[index]);}, \
        [p](){return p.size();}, \
      }; \
    }
  ACCESS()
  ACCESS(const)
  #undef ACCESS

  //! \brief The object that `this` is a proxy for in reciprocal connections.
  //! \see `Reciprocal_ptr::mine`
  Mortal_ptr<T> mine;

private:
  T* _mine() override {return mine.get();}
  const T* _mine() const override {return mine.get();}
};

}
#endif
