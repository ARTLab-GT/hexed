#include <Hanging_node_matcher.hpp>
#include <math.hpp>

namespace cartdg
{

Hanging_node_matcher::Hanging_node_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive)
: elements{fine_elements}, id{i_dim}, isp{is_positive}, n_dim{elements[0]->storage_params().n_dim},
  n_vert{elements[0]->storage_params().n_vertices()/2}
{}

void Hanging_node_matcher::match(Element::shareable_value_access access_func)
{
  // extract the indices and values of the vertices on the desired face
  int inds [4]; // for `n_dim <= 3` the required size will be <= 4
  double* value_ptrs [4];
  Eigen::VectorXd values (n_vert);
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
    int stride = custom_math::pow(2, n_dim - 1 - id);
    inds[i_vert] = i_vert/stride*stride*2 + i_vert%stride + isp*stride;
    value_ptrs[i_vert] = std::invoke(access_func, elements[i_vert]);
    values(i_vert) = value_ptrs[i_vert][inds[i_vert]];
  }
  // compute interpolation
  Eigen::MatrixXd interp_mat {
    {1. , 0. },
    {0.5, 0.5},
    {0. , 1. }
  };
  // 3[x3[x3]] array of values at corners and midpoints
  Eigen::VectorXd interpolated = custom_math::hypercube_matvec(interp_mat, values);
  // write to element vertices
  for (int i_elem = 0; i_elem < n_vert; ++i_elem) {
    // loop over vertices of this element
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
      // compute corresponding index in interpolated array
      int interp_ind = i_elem%2 + i_vert%2; // account for inner stride
      interp_ind += i_elem/2*3 + i_vert/2*3; // if 3D, account for outer stride
      value_ptrs[i_elem][inds[i_vert]] = interpolated(interp_ind);
    }
  }
}

}
