#ifndef HEXED_REFINED_FACE_HPP_
#define HEXED_REFINED_FACE_HPP_

#include <array>
#include <algorithm>

namespace hexed
{

class Refined_face
{
  public:
  double* coarse = nullptr;
  std::array<double*, 4> fine {};
  std::array<bool, 2> stretch;
  int coarse_mask = 0; //!< \note `Accessible_mesh` is supposed to set this and keep it updated
  std::array<int, 4> fine_masks {}; //!< \note `Accessible_mesh` is supposed to set this and keep it updated
  int fine_mask() {return *std::max_element(fine_masks.begin(), fine_masks.end());}
  int mask() {return std::max(coarse_mask, fine_mask());}
};

}
#endif
