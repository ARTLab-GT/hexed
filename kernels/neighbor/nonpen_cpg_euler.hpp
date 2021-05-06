#ifndef CARTDG_NONPEN_CPG_EULER_HPP_
#define CARTDG_NONPEN_CPG_EULER_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include "read_copy.hpp"
#include "write_copy.hpp"
#include "cpg_euler_hll_deformed.hpp"

namespace cartdg
{

// AUTOGENERATE
template <int n_var, int n_qpoint, int row_size>
void nonpen_cpg_euler(double* read, double* write, double* jacobian, int* i_elem,
                      int* i_dim, int* is_positive_face, int n_bc, double weight,
                      Kernel_settings& settings)
{
  const int n_face_qpoint = n_qpoint/row_size;
  const int n_dim = n_var - 2;
  double mult = settings.d_t_by_d_pos/weight;
  double heat_rat = settings.cpg_heat_rat;

  for (int i_bc = 0; i_bc < n_bc; ++i_bc)
  {
    int stride = n_face_qpoint;
    for (int i = 0; i < i_dim[i_bc]; ++i) stride /= row_size;
    bool is_positive = is_positive_face[i_bc] == 1;

    double face_r [2][n_var][n_face_qpoint];
    double face_jacobian [2][n_dim][n_dim][n_face_qpoint];
    double face_w [2][n_var][n_face_qpoint];
    read_copy<n_var, n_qpoint, row_size>(read + n_var*n_qpoint*i_elem[i_bc], &face_r[0][0][0], stride, is_positive);
    read_copy<n_dim*n_dim, n_qpoint, row_size>(jacobian + n_dim*n_dim*n_qpoint*i_elem[i_bc], &face_jacobian[0][0][0][0], stride, is_positive);

    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++ i_qpoint)
    {
      Eigen::Matrix<double, n_dim, n_dim> jac_mat;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        for (int k_dim = 0; k_dim < n_dim; ++k_dim)
        {
          jac_mat(j_dim, k_dim) = face_jacobian[0][j_dim][k_dim][i_qpoint];
        }
      }

      double normal [n_dim];
      double normal_sq = 0;
      double mom_dot_normal = 0;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        for (int k_dim = 0; k_dim < n_dim; ++k_dim)
        {
          jac_mat(k_dim, i_dim[i_bc]) = 0.;
        }
        jac_mat(j_dim, i_dim[i_bc]) = 1.;
        normal[j_dim] = jac_mat.determinant();
        mom_dot_normal += face_r[0][j_dim][i_qpoint]*normal[j_dim];
        normal_sq += normal[j_dim]*normal[j_dim];
      }
      double normal_mom [n_dim];
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        normal_mom[j_dim] = mom_dot_normal*normal[j_dim]/normal_sq;
      }

      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        face_r[1][i_var][i_qpoint] = face_r[0][i_var][i_qpoint];
      }
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        face_r[1][j_dim][i_qpoint] -= 2*normal_mom[j_dim];
        for (int k_dim = 0; k_dim < n_dim; ++k_dim)
        {
          face_jacobian[1][j_dim][k_dim][i_qpoint] = face_jacobian[0][j_dim][k_dim][i_qpoint];
        }
      }
    }

    int dims [2] {i_dim[i_bc], i_dim[i_bc]};
    bool flips [2] {!is_positive, !is_positive};
    cpg_euler_hll_deformed<n_dim, n_face_qpoint>(&face_r[0][0][0], &face_w[0][0][0], &face_jacobian[0][0][0][0], mult, dims, flips, heat_rat);
    write_copy<n_var, n_qpoint, row_size>(&face_w[0][0][0], write + n_var*n_qpoint*i_elem[i_bc], stride, is_positive);
  }
}

}

#endif
