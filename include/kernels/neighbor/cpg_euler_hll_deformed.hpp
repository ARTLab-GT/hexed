#ifndef CARTDG_CPG_EULER_HLL_DEFORMED_HPP_
#define CARTDG_CPG_EULER_HLL_DEFORMED_HPP_

#include <Eigen/Dense>

namespace cartdg
{

template<int n_dim, int n_face_qpoint>
void cpg_euler_hll_deformed(double* state_r, double* d_flux_w, double* jacobian,
                            double mult, int i_axis0, int i_axis1, double sp_heat_rat)
{
  const int n_var = n_dim + 2;
  const int face_size = n_var*n_face_qpoint;
  int i_axis_arg [] {i_axis0, i_axis1};

  for (int i_qpoint = 0; i_qpoint < n_face_qpoint; ++i_qpoint)
  {
    double flux[2][n_var];
    double wave_speed[2];
    double jacobian_det[n_face_qpoint];
    for (int i_side = 0; i_side < 2; ++i_side)
    {
      double qpoint_flux [n_dim][n_var];
      double qpoint_speed [n_dim];
      Eigen::Matrix<double, n_dim, n_dim> jac_mat;
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        #define READ(i) state_r[(i)*n_face_qpoint + i_qpoint + i_side*face_size]
        double veloc = READ(i_axis0)/READ(n_var - 2);
        double pres = 0;
        for (int j_axis = 0; j_axis < n_var - 2; ++j_axis)
        {
          qpoint_flux[i_axis][j_axis] = READ(j_axis)*veloc;
          pres += READ(j_axis)*READ(j_axis)/READ(n_var - 2);
          jac_mat(i_axis, j_axis) = jacobian[(i_axis*n_dim + j_axis)*n_face_qpoint + i_qpoint];
        }
        pres = (sp_heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
        qpoint_flux[i_axis][i_axis] += pres;
        qpoint_flux[i_axis][n_var - 2] = READ(i_axis0);
        qpoint_flux[i_axis][n_var - 1] = (READ(n_var - 1) + pres)*veloc;
        qpoint_speed[i_axis] = veloc + (2*i_side - 1)*sqrt(sp_heat_rat*pres/READ(n_var - 2));
        #undef READ
      }

      jacobian_det[i_qpoint] = jac_mat.determinant();
      for (int i_var = 0; i_var < n_var; ++i_var)
      {
        for (int i_axis = 0; i_axis < n_dim; ++i_axis)
        {
          jac_mat(i_axis, i_axis_arg[i_side]) = qpoint_flux[i_axis][i_var];
        }
        flux[i_side][i_var] = jac_mat.determinant();
      }
      for (int i_axis = 0; i_axis < n_dim; ++i_axis)
      {
        jac_mat(i_axis, i_axis_arg[i_side]) = qpoint_speed[i_axis];
      }
      wave_speed[i_side] = jac_mat.determinant();
    }

    for (int i_var = 0; i_var < n_var; ++i_var)
    {
      double num_flux;
      // Note: these conditions are different from the standard ones, but it rarely makes
      // a difference, except in the case of opposing supersonic flows, where this behavior
      // is preferable
      if      (std::min(wave_speed[0], wave_speed[1]) >= 0) num_flux = flux[0][i_var];
      else if (std::max(wave_speed[1], wave_speed[0]) <= 0) num_flux = flux[1][i_var];
      else
      {
        num_flux = (  wave_speed[1]*flux[0][i_var]
                    - wave_speed[0]*flux[1][i_var]
                    + wave_speed[0]*wave_speed[1]*(
                      state_r[i_var*n_face_qpoint + i_qpoint + face_size]
                    - state_r[i_var*n_face_qpoint + i_qpoint]))
                   / (wave_speed[1] - wave_speed[0]);
      }

      d_flux_w[i_qpoint + i_var*n_face_qpoint            ] = (flux[0][i_var] - num_flux)*mult/jacobian_det[i_qpoint];
      d_flux_w[i_qpoint + i_var*n_face_qpoint + face_size] = (num_flux - flux[1][i_var])*mult/jacobian_det[i_qpoint];
    }
  }
}

}

#endif