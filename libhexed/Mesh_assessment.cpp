#include <hexed/Mesh_assessment.hpp>

namespace hexed {

Mesh_assessment::Mesh_assessment()
: orthogonality{1}
, edge_lengths{Mat<3>::Zero()}
{}

Mesh_assessment::Mesh_assessment(next::Sequence<Mat<3>> vertices, int jac_vert, int grad_vert) {
  orthogonality = 0;
  edge_lengths.setZero();
}

}
