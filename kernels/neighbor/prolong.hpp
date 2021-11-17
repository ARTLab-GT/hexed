#ifndef CARTDG_PROLONG_HPP_
#define CARTDG_PROLONG_HPP_

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Refined_face.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void prolong(ref_face_vec& ref_faces, Basis& basis, Kernel_settings& settings)
{
  const int n_dim {n_var - 2};
  const int n_face {custom_math::(2, n_dim - 1)};

  double* coarse {ref_faces.coarse_face()};
  for (int i_face = 0; i_face < n_face; ++i_face)
  {
  }
}

}
#endif
