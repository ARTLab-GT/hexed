#ifndef HEXED_NEIGHBOR_AVG_HPP_
#define HEXED_NEIGHBOR_AVG_HPP_

#include "Vector_view.hpp"
#include "connection.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace hexed
{

template <typename element_t> bool flip(Face_connection<element_t>&);

template <> bool flip(Face_connection<Element>&)
{
  return false;
}

template <> bool flip(Face_connection<Deformed_element>& con)
{
  auto dir = con.direction();
  return dir.flip_normal(0) != dir.flip_normal(1);
}

/*
 * overwrites the data in both faces with the average of the two
 */
template <int n_dim, int row_size, typename element_t>
class Neighbor_avg : public Kernel<Face_connection<element_t>&>
{
  public:
  virtual void operator()(Sequence<Face_connection<element_t>&>& connections)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_face_qpoint = custom_math::pow(row_size, n_dim - 1);
    constexpr int face_size = n_face_qpoint*n_var;

    #pragma omp parallel for
    for (int i_con = 0; i_con < connections.size(); ++i_con)
    {
      auto& con = connections[i_con];
      double* face [2] {con.face(0), con.face(1)};
      int sign [] {1, 1 - 2*flip(con)};
      for (int i_dof = 0; i_dof < face_size; ++i_dof) {
        double avg = 0.;
        for (int i_side : {0, 1}) avg += .5*sign[i_side]*face[i_side][i_dof];
        for (int i_side : {0, 1}) face[i_side][i_dof] = sign[i_side]*avg;
      }
    }
  }
};

template <int n_dim, int row_size>
class Neighbor_avg_cartesian : public Neighbor_avg<n_dim, row_size, Element> {};

template <int n_dim, int row_size>
class Neighbor_avg_deformed : public Neighbor_avg<n_dim, row_size, Deformed_element> {};

template<>
class Kernel_traits<Neighbor_avg_cartesian>
{
  public:
  using base_t = Kernel<Face_connection<Element>&>;
};

template<>
class Kernel_traits<Neighbor_avg_deformed>
{
  public:
  using base_t = Kernel<Face_connection<Deformed_element>&>;
};

}
#endif
