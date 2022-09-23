#ifndef HEXED_NEIGHBOR_CARTESIAN_HPP_
#define HEXED_NEIGHBOR_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "hll.hpp"

namespace hexed
{

/*
 * Computes the numerical flux at the faces based on the extrapolated state data.
 * State is read from the `face`s of the face connections and the numerical flux is
 * written to the same location (the same flux is written to both faces).
 */
template <int n_dim, int row_size>
class Neighbor_cartesian : public Kernel<Face_connection<Element>&>
{
  const double heat_rat;
  public:
  Neighbor_cartesian(double heat_ratio=1.4) : heat_rat{heat_ratio} {}

  virtual void operator()(Sequence<Face_connection<Element>&>& connections)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face_qpoint = custom_math::pow(row_size, n_dim - 1);
    constexpr int face_size = n_face_qpoint*n_var;
    // all possible face normals we might need to be used in the hll flux
    double normal [n_dim][n_dim][n_face_qpoint];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint) {
          normal[i_dim][j_dim][i_qpoint] = i_dim == j_dim;
        }
      }
    }

    #pragma omp parallel for
    for (int i_con = 0; i_con < connections.size(); ++i_con)
    {
      double face [2*face_size];
      auto& con = connections[i_con];
      const int i_dim = con.direction().i_dim;
      for (int i_side : {0, 1}) {
        double* con_face = con.face(i_side);
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          face[i_side*face_size + i_face_dof] = con_face[i_face_dof];
        }
      }
      hll<n_dim, n_face_qpoint>(face, normal[i_dim][0], heat_rat);
      for (int i_side : {0, 1}) {
        double* con_face = con.face(i_side);
        for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof) {
          con_face[i_face_dof] = face[i_side*face_size + i_face_dof];
        }
      }
    }
  }
};

template<>
class Kernel_traits<Neighbor_cartesian>
{
  public:
  using base_t = Kernel<Face_connection<Element>&>;
};

}
#endif
