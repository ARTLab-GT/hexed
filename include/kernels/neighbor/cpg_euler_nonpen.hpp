#ifndef CARTDG_CPG_EULER_NONPEN_HPP_
#define CARTDG_CPG_EULER_NONPEN_HPP_

#include <Eigen/Dense>

#include "../Kernel_settings.hpp"
#include "read_copy.hpp"
#include "write_copy.hpp"

namespace cartdg
{

template <int n_var, int n_qpoint, int row_size>
void cpg_euler_nonpen(double* read, double* write, double* jacobian, int* i_elem,
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
    int sign = 2*is_positive_face[i_bc] - 1;
    bool is_positive = is_positive_face[i_bc] == 1;

    double face_r [n_var][n_face_qpoint];
    double face_jacobian[n_dim][n_dim][n_face_qpoint];
    double face_w [n_var][n_face_qpoint];
    read_copy<n_var, n_qpoint, row_size>(read + n_var*n_qpoint*i_elem[i_bc], &face_r[0][0], stride, is_positive);
    read_copy<n_dim*n_dim, n_qpoint, row_size>(jacobian + n_dim*n_dim*n_qpoint*i_elem[i_bc], &face_jacobian[0][0][0], stride, is_positive);

    for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++ i_qpoint)
    {
      Eigen::Matrix<double, n_dim, n_dim> jac_mat;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        for (int k_dim = 0; k_dim < n_dim; ++k_dim)
        {
          jac_mat(j_dim, k_dim) = face_jacobian[j_dim][k_dim][i_qpoint];
        }
      }
      double jac_det = jac_mat.determinant();
      double veloc [n_dim];
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        veloc[j_dim] = face_r[j_dim][i_qpoint]/face_r[n_var - 2][i_qpoint];
        jac_mat(j_dim, i_dim[i_bc]) = veloc[j_dim];
      }
      double veloc_normal = jac_mat.determinant()*sign;

      double kin_ener = 0.;
      for (int j_dim = 0; j_dim < n_dim; ++j_dim)
      {
        kin_ener += face_r[j_dim][i_qpoint]*veloc[j_dim];
      }
      kin_ener *= 0.5;
      double pressure = (face_r[n_var - 1][i_qpoint] - kin_ener)*(heat_rat - 1.);

      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        face_w[i_var][i_qpoint] = veloc_normal*face_r[i_var][i_qpoint]*mult/jac_det;
      }
      face_w[n_var - 1][i_qpoint] += veloc_normal*pressure*mult/jac_det;
    }

    write_copy<n_var, n_qpoint, row_size>(&face_w[0][0], write + n_var*n_qpoint*i_elem[i_bc], stride, is_positive);
  }
}

}

#endif