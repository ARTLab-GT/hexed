#ifndef CARTDG_RESTRICT_HPP_
#define CARTDG_RESTRICT_HPP_

#include <math.hpp>
#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Refined_face.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void restrict(ref_face_vec& ref_faces, Basis& basis, Kernel_settings& settings)
{
  const int n_dim {n_var - 2};
  const int n_face {custom_math::pow(2, n_dim - 1)};
  const int nfq {n_qpoint/row_size};
  Eigen::Matrix<double, row_size, row_size> restrict_mat [] {basis.restrict(0), basis.restrict(1)};

  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_ref_face = 0; i_ref_face < ref_faces[i_dim].size(); ++i_ref_face)
    {
      auto& ref_face {*ref_faces[i_dim][i_ref_face]};
      double* coarse {ref_face.coarse_face()};
      for (int i_dof = 0; i_dof < n_var*nfq; ++i_dof) coarse[i_dof] = 0.;
      for (int i_face = 0; i_face < n_face; ++i_face)
      {
        double* fine {ref_face.fine_face(i_face)};
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          double* var_face {fine + i_var*nfq};
          for (int j_dim = 0; j_dim < n_dim - 1; ++j_dim)
          {
            const int pow {n_dim - 2 - j_dim};
            const int face_stride {custom_math::pow(2, pow)};
            const int qpoint_stride {custom_math::pow(row_size, pow)};
            const int i_half {(i_face/face_stride)%2};
            for (int i_outer = 0; i_outer < nfq/(row_size*qpoint_stride); ++i_outer)
            {
              for (int i_inner = 0; i_inner < qpoint_stride; ++i_inner)
              {
                Eigen::Matrix<double, row_size, 1> row;
                for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                  row(i_qpoint) = var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner];
                }
                row = restrict_mat[i_half]*row;
                for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                  var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner] = row(i_qpoint);
                }
              }
            }
          }
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
            coarse[i_var*nfq + i_qpoint] += var_face[i_qpoint];
          }
        }
      }
    }
  }
}

}
#endif
