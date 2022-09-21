#include <connection.hpp>

namespace hexed
{

std::vector<int> face_vertex_inds(int n_dim, Con_dir<Deformed_element> direction)
{
  int n_vert = custom_math::pow(2, n_dim - 1);
  std::vector<int> inds;
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) inds.push_back(i_vert);
  // reorder as necessary
  if (direction.flip_tangential()) {
    // if there is a dimension not involved in the connection which is greater than `i_dim[0]`
    // the stride to flip is 2. otherwise it is 0.
    int stride = 1;
    if (n_dim == 3) {
      int unused_dim = 3 - direction.i_dim[0] - direction.i_dim[1];
      if (unused_dim > direction.i_dim[0]) stride = 2;
    }
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
      // this arithmetic is equivalent to swapping vertices along dimension i_dim[0]
      inds[i_vert] += stride*(1 - 2*((i_vert/stride)%2));
    }
  }
  if (direction.transpose()) std::swap(inds[1], inds[2]); // only possible for this to happen if `n_dim == 3`
  return inds;
}

std::array<std::vector<int>, 2> vertex_inds(int n_dim, Con_dir<Deformed_element> direction)
{
  // get vertices involved
  std::array<std::vector<int>, 2> inds;
  int n_vert = custom_math::pow(2, n_dim - 1);
  std::array<int, 2> strides;
  for (int i_side = 0; i_side < 2; ++i_side) {
    strides[i_side] = std::pow(2, n_dim - direction.i_dim[i_side] - 1);
    for (int i_vertex = 0; i_vertex < 2*n_vert; ++i_vertex) {
      if ((i_vertex/strides[i_side])%2 == int(direction.face_sign[i_side])) {
        inds[i_side].push_back(i_vertex);
      }
    }
  }
  std::vector<int> permutation = face_vertex_inds(n_dim, direction);
  std::vector<int> permuted;
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
    permuted.push_back(inds[1][permutation[i_vert]]);
  }
  std::swap(inds[1], permuted);
  return inds;
}

}
