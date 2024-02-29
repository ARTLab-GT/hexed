#ifndef HEXED_REFINED_FACE_HPP_
#define HEXED_REFINED_FACE_HPP_

#include <array>

namespace hexed
{

class Refined_face
{
  public:
  double* coarse = nullptr;
  std::array<double*, 4> fine {};
  std::array<bool, 2> stretch;
  bool coarse_mask = true; //!< \note `Accessible_mesh` is supposed to set this and keep it updated
  std::array<bool, 4> fine_masks {true, true, true, true}; //!< \note `Accessible_mesh` is supposed to set this and keep it updated
  bool any_fine_mask() {return std::any_of(fine_masks.begin(), fine_masks.end(), [](bool b){return b;});}
};

}
#endif
