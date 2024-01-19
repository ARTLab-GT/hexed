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
};

}
#endif
