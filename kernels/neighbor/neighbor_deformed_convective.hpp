#ifndef CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_deformed_convective(def_elem_con_vec& def_connections, Basis& basis, Kernel_settings& settings)
// "def_" (for "deformed") prepended to arguments to avoid naming conflicts in benchmark code
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  const int n_dim = n_var - 2;
  const int jac_size = n_dim*n_dim*n_face_qpoint;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;

  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < def_connections.size(); ++i_con)
  {
    double face_r [2*face_size];
    double face_jacobian [2*jac_size];
    double face_w [2*face_size];

    Deformed_elem_con connection = def_connections[i_con];
    bool is_positive [2] {connection.is_positive[0], connection.is_positive[1]};
    int i_dim [2] {connection.i_dim[0], connection.i_dim[1]};
    for (int i_side : {0, 1})
    {
      int stride = n_face_qpoint;
      for (int i = 0; i < i_dim[i_side]; ++i) stride /= row_size;
      double* stage_r = connection.element[i_side]->stage(i_read);
      double* jacobian = connection.element[i_side]->jacobian();
      read_copy<n_var, n_qpoint, row_size>(stage_r, face_r + i_side*face_size, stride, is_positive[i_side]);
      read_copy<n_dim*n_dim, n_qpoint, row_size>(jacobian, face_jacobian + i_side*jac_size, stride, is_positive[i_side]);
    }

    if ((is_positive[0] != is_positive[1]) && (i_dim[0] != i_dim[1]))
    {
      if (n_dim == 3)
      {} // FIXME
      else
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var>> face_state (face_r);
        face_state.colwise().reverseInPlace();
        Eigen::Map<Eigen::Matrix<double, row_size, n_dim*n_dim>> face_jac (face_jacobian);
        face_jac.colwise().reverseInPlace();
      }
    }

    if ((is_positive[0] != is_positive[1]) && (i_dim[0] != i_dim[1]))
    {
      if (n_dim == 3)
      {}
      else
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var>> face_state (face_w);
        face_state.colwise().reverseInPlace();
      }
    }

    for (int i_side : {0, 1})
    {
      int stride = n_face_qpoint;
      for (int i = 0; i < i_dim[i_side]; ++i) stride /= row_size;
      double* write = connection.element[i_side]->stage(i_write);
      write_copy<n_var, n_qpoint, row_size, true>(face_w + i_side*face_size, write, stride, is_positive[i_side]);
    }
  }
}

}

#endif
