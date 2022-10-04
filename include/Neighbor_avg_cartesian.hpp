#ifndef HEXED_NEIGHBOR_AVG_CARTESIAN_HPP_
#define HEXED_NEIGHBOR_AVG_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace hexed
{

/*
 * overwrites the data in both faces with the average of the two
 */
template <int n_dim, int row_size>
class Neighbor_avg_cartesian : public Kernel<Face_connection<Element>&>
{
  public:
  virtual void operator()(Sequence<Face_connection<Element>&>& connections)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face_qpoint = custom_math::pow(row_size, n_dim - 1);
    constexpr int face_size = n_face_qpoint*n_var;

    #pragma omp parallel for
    for (int i_con = 0; i_con < connections.size(); ++i_con)
    {
      auto& con = connections[i_con];
      double* face [2] {con.face(0), con.face(1)};
      for (int i_dof = 0; i_dof < face_size; ++i_dof) {
        double avg = 0.;
        for (int i_side : {0, 1}) avg += .5*face[i_side][i_dof];
        for (int i_side : {0, 1}) face[i_side][i_dof] = avg;
      }
    }
  }
};

template<>
class Kernel_traits<Neighbor_avg_cartesian>
{
  public:
  using base_t = Kernel<Face_connection<Element>&>;
};

}
#endif
