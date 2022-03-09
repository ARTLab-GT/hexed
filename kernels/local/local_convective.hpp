#ifndef CARTDG_LOCAL_CONVECTIVE_HPP_
#define CARTDG_LOCAL_CONVECTIVE_HPP_

#include <Eigen/Dense>

#include <Kernel_settings.hpp>
#include <Basis.hpp>
#include <Element.hpp>

namespace cartdg
{

// AUTOGENERATE LOOKUP BENCHMARK(cartesian, 3)
template<int n_var, int n_qpoint, int row_size>
void local_convective(elem_vec& elements, Basis& basis, Kernel_settings& settings)
{
  const Eigen::Matrix<double, row_size, row_size> diff_mat {basis.diff_mat()};
  const Eigen::Matrix<double, 2, row_size> boundary {basis.boundary()};
  Eigen::MatrixXd sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights {Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()};
  const Eigen::Matrix<double, row_size, 2> lift {inv_weights.asDiagonal()*basis.boundary().transpose()*sign};

  const double d_t_by_d_pos = settings.d_t/settings.d_pos;
  const double heat_rat = settings.cpg_heat_rat;
  const int n_face_dof = n_var*n_qpoint/row_size;
  const double rk_weight = settings.rk_weight;

  #pragma omp parallel for
  for (unsigned i_elem = 0; i_elem < elements.size(); ++i_elem)
  {
    double* state = elements[i_elem]->stage(0);
    double* rk_reference = state + n_var*n_qpoint;
    double time_rate [n_var][n_qpoint] {};
    double* face = elements[i_elem]->face();
    double* tss = elements[i_elem]->time_step_scale();

    // Compute update
    for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
         stride /= row_size, n_rows *= row_size, ++i_dim)
    {
      double* face0 = face + i_dim*2*n_face_dof;
      double* face1 = face0 + n_face_dof;
      int i_face_qpoint {0};
      for (int i_outer = 0; i_outer < n_rows; ++i_outer)
      {
        for (int i_inner = 0; i_inner < stride; ++i_inner)
        {
          // Fetch this row of data
          double row_r [n_var][row_size];
          double row_w [n_var][row_size];
          for (int i_var = 0; i_var < n_var; ++i_var) {
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              row_r[i_var][i_qpoint] = state[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
            }
          }

          // Calculate flux
          double flux [n_var][row_size];
          for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint)
          {
            #define READ(i) row_r[i][i_qpoint]
            #define FLUX(i) flux[i][i_qpoint]
            double veloc = READ(i_dim)/READ(n_var - 2);
            double pres = 0;
            for (int j_dim = 0; j_dim < n_var - 2; ++j_dim) {
              FLUX(j_dim) = READ(j_dim)*veloc;
              pres += READ(j_dim)*READ(j_dim)/READ(n_var - 2);
            }
            pres = (heat_rat - 1.)*(READ(n_var - 1) - 0.5*pres);
            FLUX(i_dim) += pres;
            FLUX(n_var - 2) = READ(i_dim);
            FLUX(n_var - 1) = (READ(n_var - 1) + pres)*veloc;
            #undef FLUX
            #undef READ
          }

          // Differentiate flux
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> f (&(flux[0][0]));
          Eigen::Map<Eigen::Matrix<double, row_size, n_var>> w (&(row_w[0][0]));
          w.noalias() = -diff_mat*f;

          // Add dimensional component to update
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            const int face_offset = i_var*n_qpoint/row_size + i_face_qpoint;
            Eigen::Matrix<double, 2, 1> boundary_values {face0[face_offset], face1[face_offset]};
            w.col(i_var).noalias() += lift*(boundary_values - boundary*f.col(i_var));
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              int offset = i_outer*stride*row_size + i_inner + i_qpoint*stride;
              time_rate[i_var][offset] += row_w[i_var][i_qpoint];
            }
          }
          ++i_face_qpoint;
        }
      }
    }

    // write the updated solution
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        const int i_dof = i_var*n_qpoint + i_qpoint;
        double updated = time_rate[i_var][i_qpoint]*d_t_by_d_pos*tss[i_qpoint] + state[i_dof];
        state[i_dof] = rk_weight*updated + (1. - rk_weight)*rk_reference[i_dof];
      }
    }
  }
}

}

#endif
