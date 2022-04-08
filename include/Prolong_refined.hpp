#ifndef CARTDG_PROLONG_REFINED_HPP_
#define CARTDG_PROLONG_REFINED_HPP_

#include "Vector_view.hpp"
#include "Refined_face.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace cartdg
{

/*
 * See `Refined_face.hpp`.
 * Transforms data from coarse face to fine faces by polynomial interpolation.
 * This interpolation is exact (and therefore conservative).
 */
template <int n_dim, int row_size>
class Prolong_refined : public Kernel<Refined_face&>
{
  const Eigen::Matrix<double, row_size, row_size> prolong_mat [2];

  public:
  Prolong_refined(const Basis& basis) : prolong_mat{basis.prolong(0), basis.prolong(1)} {}

  virtual void operator()(Sequence<Refined_face&>& ref_faces)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face = custom_math::pow(2, n_dim - 1);
    constexpr int nfq = custom_math::pow(row_size, n_dim - 1);

    #pragma omp parallel for
    for (int i_ref_face = 0; i_ref_face < ref_faces.size(); ++i_ref_face)
    {
      auto& ref_face {ref_faces[i_ref_face]};
      double* coarse {ref_face.coarse_face()};
      for (int i_face = 0; i_face < n_face; ++i_face)
      {
        double* fine {ref_face.fine_face(i_face)};
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          double* var_face {fine + i_var*nfq};
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
            var_face[i_qpoint] = coarse[i_var*nfq + i_qpoint];
          }
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
                row = prolong_mat[i_half]*row;
                for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                  var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner] = row(i_qpoint);
                }
              }
            }
          }
        }
      }
    }
  }
};

template<>
class Kernel_traits<Prolong_refined>
{
  public:
  using base_t = Kernel<Refined_face&>;
};

}
#endif