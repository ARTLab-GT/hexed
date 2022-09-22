#ifndef HEXED_DEFORMED_CARTESIAN_HPP_
#define HEXED_DEFORMED_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "hll.hpp"
#include "Face_permutation.hpp"

namespace hexed
{

/*
 * Computes the numerical flux at the faces based on the extrapolated state data.
 * State is read from the `face`s of the face connections and the numerical flux is
 * written to the same location (the same flux is written to both faces).
 */
template <int n_dim, int row_size>
class Neighbor_deformed : public Kernel<Face_connection<Deformed_element>&>
{
  const double heat_rat;
  public:
  Neighbor_deformed(double heat_ratio=1.4) : heat_rat{heat_ratio} {}

  virtual void operator()(Sequence<Face_connection<Deformed_element>&>& connections)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face_qpoint = custom_math::pow(row_size, n_dim - 1);
    constexpr int face_size = n_face_qpoint*n_var;
    #pragma omp parallel for
    for (int i_con = 0; i_con < connections.size(); ++i_con)
    {
      double face [2*face_size];
      auto& con = connections[i_con];
      double* jacobian = con.jacobian();
      double* elem_face [2];
      auto dir = con.direction();
      for (int i_side : {0, 1}) {
        elem_face[i_side] = con.face(i_side);
      }
      const int i_dim = dir.i_dim[0]; // which dimension is the face normal in face coords
      Face_permutation<n_dim, row_size> permutation(dir, face + face_size);

      // compute surface normals
      double normals [n_dim][n_face_qpoint];
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
        Eigen::Matrix<double, n_dim, n_dim> jac;
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
            jac(j_dim, k_dim) = jacobian[(j_dim*n_dim + k_dim)*n_face_qpoint + i_qpoint];
          }
        }
        auto orthonormal = custom_math::orthonormal(jac, i_dim);
        jac.col(i_dim) = orthonormal.col(i_dim);
        double jac_det = std::abs(jac.determinant()); // `abs` since surface coordinates not guaranteed to be right-hand
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          normals[j_dim][i_qpoint] = orthonormal(j_dim, i_dim)*jac_det;
        }
      }

      // fetch face state from element storage
      for (int i_side : {0, 1}) {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          face[i_side*face_size + i_face_dof] = elem_face[i_side][i_face_dof];
        }
      }

      // compute upwind flux
      permutation.match_faces();
      hll<n_dim, n_face_qpoint>(face, normals[0], heat_rat);
      permutation.restore();

      // flip normal
      for (int i_side : {0, 1}) {
        if (dir.flip_normal(i_side)) {
          for (int i_var = 0; i_var < n_var; ++i_var) {
            for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
              face[face_size*i_side + i_var*n_face_qpoint + i_qpoint] *= -1;
            }
          }
        }
      }
      // write numerical flux to element storage
      for (int i_side : {0, 1}) {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          elem_face[i_side][i_face_dof] = face[i_side*face_size + i_face_dof];
        }
      }
    }
  }
};

template<>
class Kernel_traits<Neighbor_deformed>
{
  public:
  using base_t = Kernel<Face_connection<Deformed_element>&>;
};

}
#endif
