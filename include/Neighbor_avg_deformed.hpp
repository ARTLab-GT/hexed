#ifndef HEXED_NEIGHBOR_AVG_DEFORMED_HPP_
#define HEXED_NEIGHBOR_AVG_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Face_permutation.hpp"

namespace hexed
{

/*
 * overwrites the data in both faces with the average of the two,
 * with appropriate normal-flipping and permutation
 */
template <int n_dim, int row_size>
class Neighbor_avg_deformed : public Kernel<Face_connection<Deformed_element>&>
{
  bool f;
  int side;

  public:
  Neighbor_avg_deformed(bool flip, int which_side = 0) : f{flip}, side{which_side} {}

  virtual void operator()(Sequence<Face_connection<Deformed_element>&>& connections)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face_qpoint = custom_math::pow(row_size, n_dim - 1);
    constexpr int face_size = n_face_qpoint*n_var;

    #pragma omp parallel for
    for (int i_con = 0; i_con < connections.size(); ++i_con)
    {
      auto& con = connections[i_con];
      double* face [2] {con.face(0), con.face(1)};
      auto dir = con.direction();
      int sign [] {1, 1 - 2*(f && (dir.flip_normal(0) != dir.flip_normal(1)))};
      Face_permutation<n_dim, row_size> permutation(dir, face[1]);
      permutation.match_faces();
      for (int i_dof = 0; i_dof < face_size; ++i_dof) {
        double avg = sign[side]*face[side][i_dof];
        for (int i_side : {0, 1}) face[i_side][i_dof] = sign[i_side]*avg;
      }
      permutation.restore();
    }
  }
};

}
#endif
