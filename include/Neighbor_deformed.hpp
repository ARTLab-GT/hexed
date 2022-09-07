#ifndef CARTDG_DEFORMED_CARTESIAN_HPP_
#define CARTDG_DEFORMED_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "hll.hpp"
#include "Surface_rotation.hpp"
#include "Face_permutation.hpp"

namespace cartdg
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
      Surface_rotation<n_dim, row_size> rotation {jacobian, i_dim};
      Face_permutation<n_dim, row_size> permutation(dir, face + face_size);

      // fetch face state from element storage
      for (int i_side : {0, 1}) {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          face[i_side*face_size + i_face_dof] = elem_face[i_side][i_face_dof];
        }
      }
      permutation.match_faces();
      // rotate momentum into surface coordinates
      for (int i_side : {0, 1}) {
        rotation.to_surface(face + i_side*face_size);
      }

      // compute upwind flux
      hll<n_dim, n_face_qpoint>(face, i_dim, heat_rat);

      // rotate momentum back into physical coordinates
      for (int i_side : {0, 1}) {
        rotation.from_surface(face + i_side*face_size);
      }
      // multiply flux by Jacobian determinant. `2*i_var` to cover both sides
      for (int i_var = 0; i_var < 2*n_var; ++i_var) {
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
          face[i_var*n_face_qpoint + i_qpoint] *= rotation.jacobian_determinant(i_qpoint);
        }
      }
      permutation.restore();
      // re-flip normal
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
