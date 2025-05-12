#ifndef HEXED_MESH_ASSESSMENT_HPP_
#define HEXED_MESH_ASSESSMENT_HPP_

#include "math.hpp"
#include "Sequence.hpp"

namespace hexed {

struct Mesh_assessment {
  Mesh_assessment();
  Mesh_assessment(next::Sequence<Mat<3>> vertices, int jac_vert, int grad_vert);
  Mat<3> orthogonality;
  Mat<3> edge_lengths;
  Mat<3, 3> grad_orth;
  Mat<3, 3> grad_lengths;
};

}
#endif
