#ifndef HEXED_FACE_HPP_
#define HEXED_FACE_HPP_

#include "reciprocal.hpp"
#include "Storage_params.hpp"
#include "Neighbor_connection.hpp"

namespace hexed {

class Element;

class Face : public Mortal {
  public:
  Face(Storage_params, int i_dim, int sign);
  inline int i_dim() const {return _i_dim;}
  inline int sign() const {return _sign;}
  void associate(Element&);
  inline Element* element() {return _element.get();}
  inline bool associated() const {return _element;}
  void connect(Reciprocal_ptr<Neighbor_connection, Face>&);
  inline Neighbor_connection* neighbor_connection() {return _neighbor_connection.get();}
  inline bool connected() const {return _neighbor_connection;}

  private:
  int _i_dim;
  int _sign;
  Mortal_ptr<Element> _element;
  Reciprocal_ptr<Face, Neighbor_connection> _neighbor_connection;
};

}
#include "Element.hpp"
#endif
