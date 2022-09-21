#ifndef HEXED_LOCAL_DEFORMED_HPP_
#define HEXED_LOCAL_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "Deformed_element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"

namespace hexed
{

/*
 * Computes the local update for one Runge-Kutta stage.
 * See `Local_cartesian_dynamic`.
 */
template <int n_dim, int row_size>
class Local_deformed : public Kernel<Deformed_element&>
{
  double dt;
  double rkw;
  const double heat_rat;

  const Eigen::Matrix<double, row_size, row_size> diff_mat;
  const Eigen::Matrix<double, 2, row_size> boundary;
  const Eigen::MatrixXd sign {{1, 0}, {0, -1}};
  const Eigen::Matrix<double, row_size, 1> inv_weights;
  const Eigen::Matrix<double, row_size, 2> lift;

  public:
  Local_deformed(const Basis& basis, double d_time, double rk_weight, double heat_ratio=1.4) :
    dt{d_time},
    rkw{rk_weight},
    heat_rat{heat_ratio},
    diff_mat{basis.diff_mat()},
    boundary{basis.boundary()},
    inv_weights{Eigen::Array<double, row_size, 1>::Constant(1.)/basis.node_weights().array()},
    lift{inv_weights.asDiagonal()*basis.boundary().transpose()*sign}
  {}

  virtual void operator()(Sequence<Deformed_element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      auto& elem {elements[i_elem]};
      double* state  = elem.stage(0);
      double* rk_reference = state + n_var*n_qpoint;
      double time_rate [n_var][n_qpoint] {};
      double* jacobian = elem.jacobian();
      double* face = elem.face();
      double* tss = elem.time_step_scale();
      double flux [n_dim][n_var][n_qpoint];
      const double d_t_by_d_pos = dt/elem.nominal_size();

      // precompute flux
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          #define READ(i) state[(i)*n_qpoint + i_qpoint]
          #define FLUX(i) flux[i_dim][i][i_qpoint]
          double veloc = READ(i_dim)/READ(n_var - 2);
          double pres = 0;
          for (int j_dim = 0; j_dim < n_var - 2; ++j_dim)
          {
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
      }

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
            // compute flux derivative
            Eigen::Matrix<double, row_size, n_dim*n_var> row_f;
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              for (int i_var = 0; i_var < n_var; ++i_var) {
                for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                  row_f(i_qpoint, j_dim*n_var + i_var) = flux[j_dim][i_var][i_outer*stride*row_size + i_inner + i_qpoint*stride];
                }
              }
            }
            Eigen::Matrix<double, row_size, n_dim*n_var> d_flux {-diff_mat*row_f};

            // fetch jacobian
            Eigen::Matrix<double, row_size, n_dim*n_dim> row_jac; // Jacobian entries in column-major order
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
                for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                  // fetch from row-major element storage
                  row_jac(i_qpoint, j_dim + n_dim*k_dim) = jacobian[(j_dim*n_dim + k_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_qpoint*stride];
                }
              }
            }
            // transform flux to reference space & compute Jacobian det
            Eigen::Matrix<double, row_size, n_var> row_w;
            Eigen::Array<double, row_size, 1> jac_det;
            for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
              // fetch qpoint Jacobian
              Eigen::Matrix<double, n_dim, n_dim> jac;
              for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) jac(i_jac) = row_jac(i_qpoint, i_jac);
              // compute determinant
              jac_det(i_qpoint) = jac.determinant();
              // compute reference flux derivative
              for (int i_var = 0; i_var < n_var; ++i_var) {
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  jac(j_dim, i_dim) = d_flux(i_qpoint, j_dim*n_var + i_var);
                }
                row_w(i_qpoint, i_var) = jac.determinant();
              }
            }

            // fetch numerical flux at faces
            Eigen::Matrix<double, n_var, 2> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              const int face_offset = i_var*n_qpoint/row_size + i_face_qpoint;
              boundary_values(i_var, 0) = face0[face_offset];
              boundary_values(i_var, 1) = face1[face_offset];
            }
            // extrapolate (reference) flux to boundaries
            Eigen::Matrix<double, 2, n_dim*n_dim> extrap_jac {boundary*row_jac};
            Eigen::Matrix<double, 2, n_dim*n_var> extrap_flux {boundary*row_f};
            for (int i_side : {0, 1}) {
              Eigen::Matrix<double, n_dim, n_dim> jac;
              for (int i_jac = 0; i_jac < n_dim*n_dim; ++i_jac) jac(i_jac) = extrap_jac(i_side, i_jac);
              Eigen::Matrix<double, n_var, 1> normal_flux;
              for (int i_var = 0; i_var < n_var; ++i_var) {
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  jac(j_dim, i_dim) = extrap_flux(i_side, j_dim*n_var + i_var);
                }
                normal_flux(i_var) = jac.determinant();
              }
              boundary_values.col(i_side) -= normal_flux;
            }

            // Add dimensional component to update
            for (int i_var = 0; i_var < n_var; ++i_var) {
              row_w.col(i_var).noalias() += lift*boundary_values.row(i_var).transpose();
              row_w.col(i_var).array() /= jac_det;
              for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                int offset = i_outer*stride*row_size + i_inner + i_qpoint*stride;
                time_rate[i_var][offset] += row_w(i_qpoint, i_var);
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
          state[i_dof] = rkw*updated + (1. - rkw)*rk_reference[i_dof];
        }
      }
    }
  }
};

template<>
class Kernel_traits<Local_deformed>
{
  public:
  using base_t = Kernel<Deformed_element&>;
};

}
#endif
