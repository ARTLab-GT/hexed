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
  protected:
  virtual void compute_flux(double*, double*) = 0;

  public:
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
      double face_nrml [n_dim*n_face_qpoint];
      double* elem_face [2];
      auto dir = con.direction();
      for (int i_side : {0, 1}) {
        elem_face[i_side] = con.face(i_side);
      }
      Face_permutation<n_dim, row_size> permutation(dir, face + face_size);

      // fetch face state from element storage
      for (int i_side : {0, 1}) {
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          face[i_side*face_size + i_face_dof] = elem_face[i_side][i_face_dof];
        }
      }
      // flip normal if necessary
      {
        int sign = 1 - 2*dir.flip_normal(0);
        double* n = con.normal(0);
        for (int i_normal = 0; i_normal < n_dim*n_face_qpoint; ++i_normal) {
          face_nrml[i_normal] = sign*n[i_normal];
        }
      }

      // compute upwind flux
      permutation.match_faces();
      compute_flux(face, face_nrml);
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

template <int n_dim, int row_size>
class Neighbor_inviscid_deformed : public Neighbor_deformed<n_dim, row_size>
{
  protected:
  const double heat_rat;
  virtual void compute_flux(double* face, double* face_nrml)
  {
    hll::inviscid<n_dim, custom_math::pow(row_size, n_dim - 1)>(face, face_nrml, heat_rat);
  }

  public:
  Neighbor_inviscid_deformed(double heat_ratio=1.4) : heat_rat{heat_ratio} {}
};

template <int n_dim, int row_size>
class Neighbor_advection_deformed : public Neighbor_deformed<n_dim, row_size>
{
  protected:
  virtual void compute_flux(double* face, double* face_nrml)
  {
    hll::advection<n_dim, custom_math::pow(row_size, n_dim - 1)>(face, face_nrml);
  }
};

}
#endif
