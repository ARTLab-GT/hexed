#include <hexed/Mesh_assessment.hpp>

namespace hexed {

Mesh_assessment::Mesh_assessment()
: orthogonality{Mat<3>::Ones()}
, edge_lengths{Mat<3>::Zero()}
, grad_orth{Mat<3, 3>::Zero()}
, grad_lengths{Mat<3, 3>::Zero()}
{}

Mesh_assessment::Mesh_assessment(next::Sequence<Mat<3>> vertices, int jac_vert, int grad_vert) {
  int n_dim = math::log(2, vertices.size());
  HEXED_ASSERT(math::pow(2, n_dim) == vertices.size(), "Size of `vertices` is not a power of 2.");
  edge_lengths.setZero();
  Mat<3, 3> edges = Mat<3, 3>::Identity();
  int edge_sensitivity [3] {};
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    int stride = math::pow(2, n_dim - 1 - i_dim);
    int start = jac_vert - jac_vert/stride%2*stride;
    int end = start + stride;
    Mat<3> edge = vertices[end] - vertices[start];
    edge_lengths(i_dim) = edge.norm();
    edges(all, i_dim) = edge/edge_lengths(i_dim);
    if (grad_vert == start) edge_sensitivity[i_dim] = -1;
    if (grad_vert ==   end) edge_sensitivity[i_dim] =  1;
  }
  orthogonality.setOnes();
  orthogonality(0) = edges.determinant();
  grad_orth.setZero();
  grad_lengths.setZero();
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    Mat<3> plane_normal = edges(all, (i_dim + 1)%3).cross(edges(all, (i_dim + 2)%3)).normalized();
    orthogonality(i_dim) = edges(all, i_dim).dot(plane_normal);
    Mat<3> dim_grad = edges(all, (i_dim + 1)%3).cross(edges(all, (i_dim + 2)%3));
    dim_grad -= dim_grad.dot(edges(all, i_dim))*edges(all, i_dim);
    grad_orth(all, 0) += dim_grad*edge_sensitivity[i_dim]/edge_lengths(i_dim);
    grad_lengths(all, i_dim) = edge_sensitivity[i_dim]*edges(all, i_dim).transpose();
  }
}

}
