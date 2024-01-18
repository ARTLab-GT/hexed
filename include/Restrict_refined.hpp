#ifndef HEXED_RESTRICT_REFINED_HPP_
#define HEXED_RESTRICT_REFINED_HPP_

#include "Sequences.hpp"
#include "connection.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace hexed
{

/*
 * See `Refined_face.hpp`.
 * Transforms data from fine faces to coarse face by orthogonal projection.
 * This projection is conservative (in the sense that the integral of the face data does not change)
 * but of course not always exact.
 */
template <int n_dim, int row_size>
class Restrict_refined : public Kernel<Refined_face&>
{
  const Eigen::Matrix<double, row_size, row_size> restrict_mat [2];
  bool scl;
  bool off;

  public:
  Restrict_refined(const Basis& basis, bool scale = true, bool offset = false) :
    restrict_mat{basis.restrict(0), basis.restrict(1)},
    scl{scale},
    off{offset}
  {}

  virtual void operator()(Sequence<Refined_face&>& ref_faces)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face = math::pow(2, n_dim - 1);
    constexpr int nfq = math::pow(row_size, n_dim - 1);

    #pragma omp parallel for
    for (int i_ref_face = 0; i_ref_face < ref_faces.size(); ++i_ref_face)
    {
      auto& ref_face {ref_faces[i_ref_face]};
      double* coarse {ref_face.coarse + 2*off*(n_dim + 2)*nfq};
      for (int i_dof = 0; i_dof < n_var*nfq; ++i_dof) coarse[i_dof] = 0.;
      auto str = ref_face.stretch;
      // update number of faces to reflect any face stretching
      int nf = n_face;
      for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim) nf /= 1 + str[i_dim];
      for (int i_face = 0; i_face < nf; ++i_face)
      {
        double* fine {ref_face.fine[i_face] + 2*off*(n_dim + 2)*nfq};
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          double* var_face {fine + i_var*nfq};
          for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim)
          {
            if (str[i_dim])
            {
              for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) var_face[i_qpoint] /= 1 + scl;
            }
            else
            {
              const int pow {n_dim - 2 - i_dim};
              const int face_stride {str[n_dim - 2] ? 1 : math::pow(2, pow)};
              const int qpoint_stride {math::pow(row_size, pow)};
              const int i_half {(i_face/face_stride)%2}; // is this face covering the upper or lower half of the coarse face with respect to the current dimension?
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
          }
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
            coarse[i_var*nfq + i_qpoint] += var_face[i_qpoint];
          }
        }
      }
    }
  }
};

}
#endif
