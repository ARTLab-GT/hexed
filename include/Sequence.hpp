#ifndef HEXED_SEQUENCE_HPP_
#define HEXED_SEQUENCE_HPP_

namespace hexed
{

/*! \brief An interface for general sequence-type containers which supports access to elements
 * but intentionally doesn't support insertion or removal of elements.
 * \details This is useful for classes providing limited public access to otherwise private variables.
 */
template<typename T>
class Sequence
{
  public:
  virtual int size() = 0;
  virtual T operator[](int index) = 0;
};

}
#endif
