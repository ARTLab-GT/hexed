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
  double mult = settings.d_t_by_d_pos/basis.node_weights()(0);
  double heat_rat = settings.cpg_heat_rat;
  const int i_read = settings.i_read;
  const int i_write = settings.i_write;

  const int face_size = n_face_qpoint*n_var;
  for (unsigned stride = n_face_qpoint, i_dim = 0; stride > 0; stride /= row_size, ++i_dim)
  {
    #pragma omp parallel for
    for (unsigned i_con = 0; i_con < def_reg_cons[i_dim].size(); ++i_con)
    {
      double face_r [2*face_size];
      double face_w [2*face_size];
      Deformed_element* def = def_reg_cons[i_dim][i_con].first;
      Element*          reg = def_reg_cons[i_dim][i_con].second;

      read_copy<n_var, n_qpoint, row_size>(def->stage(i_read), face_r            , stride, 1);
      read_copy<n_var, n_qpoint, row_size>(reg->stage(i_read), face_r + face_size, stride, 0);

      hll_cpg_euler<n_var - 2, n_face_qpoint>(face_r, face_w, mult, i_dim, heat_rat);

      write_copy<n_var, n_qpoint, row_size>(face_w            , def->stage(i_write), stride, 1);
      write_copy<n_var, n_qpoint, row_size>(face_w + face_size, reg->stage(i_write), stride, 0);
    }
  }
}

}
#endif
