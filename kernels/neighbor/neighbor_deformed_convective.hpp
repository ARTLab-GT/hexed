#ifndef CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_DEFORMED_CONVECTIVE_HPP_

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(deformed, 3)
template<int n_var, int n_qpoint, int row_size>
void neighbor_deformed_convective(def_elem_con_vec& def_connections, Basis& basis, Kernel_settings& settings)
// "def_" (for "deformed") prepended to arguments to avoid naming conflicts in benchmark code
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int face_size = n_face_qpoint*n_var;
  double heat_rat = settings.cpg_heat_rat;

  #pragma omp parallel for
  for (unsigned i_con = 0; i_con < def_connections.size(); ++i_con)
  {
    double face_r [2*face_size];
    double face_w [2*face_size];
    Deformed_elem_con con = def_connections[i_con];

    for (int i_side : {0, 1})
    {
      Face_index face_ind = con.face_index(i_side);
      double* read = face_ind.element->face() + (face_ind.i_dim*2 + face_ind.is_positive)*face_size;
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
      {
        face_r[i_side*face_size + i_face_dof] = read[i_face_dof];
      }
    }

    bool same_is_positive = (con.face_index(0).is_positive == con.face_index(1).is_positive);
    bool same_dim = (con.face_index(0).i_dim == con.face_index(1).i_dim);
    bool expected_positive [] {1, 0};
    for (int i_side : {0, 1})
    {
      if (con.face_index(i_side).is_positive != expected_positive[i_side])
      {
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
        {
          face_r[face_size*i_side + con.face_index(i_side).i_dim*n_face_qpoint + i_qpoint] *= -1;
        }
      }
    }
    // FIXME: make this work in 3D
    if (!same_dim)
    {
      if (!same_is_positive)
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_face_qpoint/row_size>> rows {face_r + face_size};
        rows.colwise().reverseInPlace();
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
        {
          face_r[face_size + con.face_index(0).i_dim*n_face_qpoint + i_qpoint] *= -1;
        }
      }
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
      {
        std::swap(face_r[face_size + con.face_index(0).i_dim*n_face_qpoint + i_qpoint],
                  face_r[face_size + con.face_index(1).i_dim*n_face_qpoint + i_qpoint]);
      }
    }

    hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, 1., con.face_index(0).i_dim, heat_rat);

    if (!same_dim)
    {
      for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
      {
        std::swap(face_w[face_size + con.face_index(0).i_dim*n_face_qpoint + i_qpoint],
                  face_w[face_size + con.face_index(1).i_dim*n_face_qpoint + i_qpoint]);
      }
      if (!same_is_positive)
      {
        Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_face_qpoint/row_size>> rows {face_w + face_size};
        rows.colwise().reverseInPlace();
        for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
        {
          face_w[face_size + con.face_index(0).i_dim*n_face_qpoint + i_qpoint] *= -1;
        }
      }
    }
    for (int i_side : {0, 1})
    {
      Face_index face_ind = con.face_index(i_side);
      if (face_ind.is_positive != expected_positive[i_side])
      {
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          if (i_var != face_ind.i_dim)
          {
            for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
            {
              face_w[face_size*i_side + i_var*n_face_qpoint + i_qpoint] *= -1;
            }
          }
        }
      }
    }

    for (int i_side : {0, 1})
    {
      Face_index face_ind = con.face_index(i_side);
      double* write = face_ind.element->face() + (face_ind.i_dim*2 + face_ind.is_positive)*face_size;
      for (int i_face_dof = 0; i_face_dof < face_size; ++i_face_dof)
      {
        write[i_face_dof] = face_w[i_side*face_size + i_face_dof];
      }
    }
  }
}

}

#endif
