#ifndef HEXED_FACE_HPP_
#define HEXED_FACE_HPP_

#include "reciprocal.hpp"
#include "Storage_params.hpp"

namespace hexed {

class Element;

class Face : public Mortal {
  public:
  Face(Storage_params, int i_dim, int sign);
  inline int i_dim() const {return _i_dim;}
  inline int sign() const {return _sign;}
  void associate(Element&);
  inline Element* element() {return _element.get();}

  private:
  int _i_dim;
  int _sign;
  Mortal_ptr<Element> _element;
};

}
#include "Element.hpp"
#endif
