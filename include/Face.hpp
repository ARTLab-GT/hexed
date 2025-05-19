#ifndef HEXED_FACE_HPP_
#define HEXED_FACE_HPP_

#include "reciprocal.hpp"
#include "Storage_params.hpp"
#include "Neighbor_connection.hpp"

namespace hexed {

class Element;
class Face_refinement;

class Face : public Mortal {
  public:
  Face(Storage_params, int i_dim, int sign);
  inline int i_dim() const {return _i_dim;}
  inline int sign() const {return _sign;}
  inline Storage_params storage_params() const {return _params;}
  void associate(Element&);
  void associate(Face_refinement&);
  inline Element* element() {return _element.get();}
  inline Face_refinement* face_ref_coarse() {return _face_ref_coarse.get();}
  inline bool associated() const {return _element || _face_ref_coarse;}
  void connect(Reciprocal_ptr<Neighbor_connection, Face>&);
  void connect(Reciprocal_ptr<Face_refinement, Face>&);
  inline Neighbor_connection* neighbor_connection() {return _neighbor_connection.get();}
  inline Face_refinement* face_ref_fine() {return _face_ref_fine.get();}
  inline bool connected() const {return _neighbor_connection || _face_ref_fine;}

  private:
  Storage_params _params;
  int _i_dim;
  int _sign;
  Mortal_ptr<Element> _element;
  Mortal_ptr<Face_refinement> _face_ref_coarse;
  Reciprocal_ptr<Face, Neighbor_connection> _neighbor_connection;
  Reciprocal_ptr<Face, Face_refinement> _face_ref_fine;
};

}
#include "Element.hpp"
#include "Face_refinement.hpp"
#endif
