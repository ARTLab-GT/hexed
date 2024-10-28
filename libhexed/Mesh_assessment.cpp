#include <hexed/Mesh_assessment.hpp>

namespace hexed {

Mesh_assessment::Mesh_assessment()
: orthogonality{1}
, edge_lengths{Mat<3>::Zero()}
{}

Mesh_assessment::Mesh_assessment(next::Sequence<Mat<3>> vertices, int jac_vert, int grad_vert) {
  int n_dim = math::log(2, vertices.size());
  HEXED_ASSERT(math::pow(2, n_dim) == vertices.size(), "Size of `vertices` is not a power of 2.");
  edge_lengths.setZero();
  Mat<3, 3> normalized_edges = Mat<3, 3>::Identity();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    int stride = math::pow(2, n_dim - 1 - i_dim);
    int start = jac_vert - jac_vert/stride%2*stride;
    int end = start + stride;
    Mat<3> edge = vertices[end] - vertices[start];
    edge_lengths(i_dim) = edge.norm();
    normalized_edges(all, i_dim) = edge/edge_lengths(i_dim);
  }
  orthogonality = normalized_edges.determinant();
}

}
