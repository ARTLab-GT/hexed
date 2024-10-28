#ifndef HEXED_MESH_ASSESSMENT_HPP_
#define HEXED_MESH_ASSESSMENT_HPP_

#include "math.hpp"
#include "Sequence.hpp"

namespace hexed {

struct Mesh_assessment {
  Mesh_assessment();
  Mesh_assessment(next::Sequence<Mat<3>> vertices, int jac_vert, int grad_vert);
  double orthogonality;
  Mat<3> edge_lengths;
};

}
#endif
