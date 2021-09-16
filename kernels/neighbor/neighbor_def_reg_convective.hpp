#ifndef CARTDG_NEIGHBOR_DEF_REG_CONVECTIVE_HPP_
#define CARTDG_NEIGHBOR_DEF_REG_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Deformed_element.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "hll_cpg_euler.hpp"

namespace cartdg
{

// AUTOGENERATE LOOKUP
template<int n_var, int n_qpoint, int row_size>
void neighbor_def_reg_convective(def_reg_con_vec& def_reg_cons, Basis& basis, Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int n_dim = n_var - 2;
  double mult = settings.d_t_by_d_pos/basis.node_weights()(0);
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;

  const int face_size = n_face_qpoint*n_var;
  for (int positive : {0, 1})
  {
    for (unsigned stride = n_face_qpoint, i_dim = 0; stride > 0; stride /= row_size, ++i_dim)
    {
      const int i_jac = i_dim*(n_dim + 1)*n_qpoint;
      const int i_vec = positive*n_dim + i_dim;
      #pragma omp parallel for
      for (unsigned i_con = 0; i_con < def_reg_cons[i_vec].size(); ++i_con)
      {
        double face_r [2*face_size];
        double face_w [2*face_size];
        double jacobian [n_face_qpoint];
        Deformed_element* def = def_reg_cons[i_vec][i_con].first;
        Element*          reg = def_reg_cons[i_vec][i_con].second;
        Element* elements [2] {positive?def:reg, positive?reg:def};

        read_copy<n_var, n_qpoint, row_size>(elements[0]->stage(i_read), face_r            , stride, 1);
        read_copy<n_var, n_qpoint, row_size>(elements[1]->stage(i_read), face_r + face_size, stride, 0);
        read_copy<1, n_qpoint, row_size>(def->jacobian() + i_jac, jacobian, stride, 1);

        hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, mult, i_dim, heat_rat);
        for (int i_var = 0; i_var < n_var; ++i_var)
        {
          int var_start = (i_var + (1 - positive)*n_var)*n_face_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
          {
            face_w[var_start + i_qpoint] /= jacobian[i_qpoint];
          }
        }

        write_copy<n_var, n_qpoint, row_size>(face_w            , elements[0]->stage(i_write), stride, 1);
        write_copy<n_var, n_qpoint, row_size>(face_w + face_size, elements[1]->stage(i_write), stride, 0);
      }
    }
  }
}

}
#endif
