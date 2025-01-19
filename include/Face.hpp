#ifndef HEXED_FACE_HPP_
#define HEXED_FACE_HPP_

#include "reciprocal.hpp"
#include "Element.hpp"

namespace hexed {

class Face : public Mortal {
  public:
  Face(Storage_params, int i_dim, int sign);
  inline int i_dim() const {return _i_dim;}
  inline int sign() const {return _sign;}
  void associate(Element&);
  void dissociate();
  inline Element* element() {return _element.get();}

  private:
  int _i_dim;
  int _sign;
  Reciprocal_ptr<Face, Element> _element;
};

}
#endif
