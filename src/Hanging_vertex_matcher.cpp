#include <Hanging_vertex_matcher.hpp>
#include <math.hpp>

namespace cartdg
{

Hanging_vertex_matcher::Hanging_vertex_matcher(std::vector<Element*> fine_elements, int i_dim, bool is_positive, std::array<bool, 2> stretch)
: elements{fine_elements}, id{i_dim}, isp{is_positive}, n_dim{elements[0]->storage_params().n_dim},
  n_vert{elements[0]->storage_params().n_vertices()/2},
  str{stretch}
{}

void Hanging_vertex_matcher::match(Element::vertex_value_access access_func)
{
  // extract the indices and values of the vertices on the desired face
  int inds [4]; // for `n_dim <= 3` the required size will be <= 4
  Eigen::VectorXd values (n_vert);
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
    int stride = custom_math::pow(2, n_dim - 1 - id);
    inds[i_vert] = i_vert/stride*stride*2 + i_vert%stride + isp*stride;
    values(i_vert) = std::invoke(access_func, elements[custom_math::stretched_ind(n_dim, i_vert, str)], inds[i_vert]);
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
  for (int i_elem = 0; i_elem < int(elements.size()); ++i_elem) {
    // loop over vertices of this element
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
      // compute corresponding index in interpolated array
      int interp_ind = (i_elem*!str[n_dim - 2])%2 + (i_vert%2)*(1 + str[n_dim - 2]); // account for inner stride
      if (n_dim == 3) interp_ind += ((i_elem*!str[0])/(1 + !str[1]) + i_vert/2*(1 + str[0]))*3; // account for outer stride
      std::invoke(access_func, elements[i_elem], inds[i_vert]) = interpolated(interp_ind);
    }
  }
}

}
