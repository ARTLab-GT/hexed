#ifndef HEXED_FACE_REFINEMENT_HPP_
#define HEXED_FACE_REFINEMENT_HPP_

#include "reciprocal.hpp"
#include "Face.hpp"

namespace hexed {

class Face_refinement : public Mortal {
  public:
  Face_refinement(Face&, int split_dim);
  inline Face& coarse() {return _coarse.value();}
  inline const Face& coarse() const {return _coarse.value();}
  inline bool is_deformed() {return _fine0.is_deformed();}
  int split_dim() const {return _split_dim;}
  inline std::array<Face*, 2> fine() {return {&_fine0, &_fine1};}
  inline std::array<const Face*, 2> fine() const {return {&_fine0, &_fine1};}
  inline bool alive() const {return _coarse;}
  std::array<std::vector<Element*>, 2> elements();
  Connection_direction get_direction();
  Kernel_face_refinement kernel_face_refinement();

  private:
  std::pair<Connection_direction, bool> _dir_reverse();
  Reciprocal_ptr<Face_refinement, Face> _coarse;
  int _split_dim;
  Face _fine0;
  Face _fine1;
};

}
#endif
