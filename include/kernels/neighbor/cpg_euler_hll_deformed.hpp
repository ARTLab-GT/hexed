#ifndef CARTDG_CPG_EULER_HLL_DEFORMED_HPP_
#define CARTDG_CPG_EULER_HLL_DEFORMED_HPP_

#include <Eigen/Dense>

namespace cartdg
{

template<int n_dim, int n_face_qpoint>
void cpg_euler_hll_deformed(double* state_r, double* d_flux_w, double* jacobian,
                            double mult, int i_axis_arg [2], bool flip [2], double sp_heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[2][n_var];
    double jacobian_det[2];

    double velocity[2];
    double sound_speed[2];
    for (int i_side = 0; i_side < 2; ++i_side)
    {
      #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
      double pres = 0;
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        pres += READ(i_axis)*READ(i_axis)/READ(n_var - 2);
      }
      pres = (sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
      sound_speed[i_side] = std::sqrt(sp_heat_rat*pres/READ(n_var - 2));

      double qpoint_flux [n_dim][n_var];
      double veloc [n_dim];
      Eigen::Matrix<double, n_dim, n_dim> jac_mat;
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        veloc[i_axis] = READ(i_axis)/READ(n_var - 2);
        for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
        {
          qpoint_flux[i_axis][j_axis] = READ(j_axis)*veloc[i_axis];
          jac_mat(i_axis, j_axis) = jacobian[(i_axis*n_dim + j_axis
                                              + i_side*n_dim*n_dim)*n_face_qpoint + i_qpoint];
        }
        qpoint_flux[i_axis][i_axis] += pres;
        qpoint_flux[i_axis][n_var - 2] = READ(i_axis);
        qpoint_flux[i_axis][n_var - 1] = (READ(n_var - 1) + pres)*veloc[i_axis];
      }
      #undef READ

      const int normal_dir = flip[i_side] ? -1 : 1;
      jacobian_det[i_side] = jac_mat.determinant();
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        for (int i_axis = 0; i_axis < n_dim; ++i_axis)
        {
          jac_mat(i_axis, i_axis_arg[i_side]) = qpoint_flux[i_axis][i_var];
        }
        flux[i_side][i_var] = jac_mat.determinant()*normal_dir;
      }
      double normal_magnitude = 0.;
      jac_mat.col(i_axis_arg[i_side]) = Eigen::Matrix<double, n_dim, 1>::Zero();
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        jac_mat(i_axis, i_axis_arg[i_side]) = 1;
        const double det = jac_mat.determinant();
        normal_magnitude += det*det;
        jac_mat(i_axis, i_axis_arg[i_side]) = 0;
      }
      normal_magnitude = std::sqrt(normal_magnitude);
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        jac_mat(i_axis, i_axis_arg[i_side]) = veloc[i_axis];
      }
      veloc[i_side] = jac_mat.determinant()*normal_dir;
    }

    double wave_speed [2];
    wave_speed[0] = std::min(velocity[0] - sound_speed[0], velocity[1] - sound_speed[1]);
    wave_speed[1] = std::max(velocity[0] + sound_speed[0], velocity[1] + sound_speed[1]);
    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      const int i = i_var*n_face_qpoint + i_qpoint;
      double num_flux;
      if      (wave_speed[0] >= 0) num_flux = flux[0][i_var];
      else if (wave_speed[1] <= 0) num_flux = flux[1][i_var];
      else
      {
        num_flux = (wave_speed[1]*flux[0][i_var] - wave_speed[0]*flux[1][i_var]
                    + wave_speed[0]*wave_speed[1]*(state_r[i + face_size] - state_r[i]))
                   / (wave_speed[1] - wave_speed[0]);
      }
      d_flux_w[i            ] = (flux[0][i_var] - num_flux)*mult/jacobian_det[0];
      d_flux_w[i + face_size] = (num_flux - flux[1][i_var])*mult/jacobian_det[1];
    }
  }
}

}

#endif