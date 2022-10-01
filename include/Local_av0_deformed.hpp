#ifndef HEXED_LOCAL_AV0_DEFORMED_HPP_
#define HEXED_LOCAL_AV0_DEFORMED_HPP_

#include "Vector_view.hpp"
#include "Deformed_element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"

namespace hexed
{

/*
 * Computes the local update for the first phase of the artificial viscosity operator (for deformed elements).
 * Requires that faces contain the shared LDG state.
 */
template <int n_dim, int row_size>
class Local_av0_deformed : public Kernel<Deformed_element&>
{
  Derivative<row_size> derivative;
  double dt;
  double rkw;

  public:
  Local_av0_deformed(const Basis& basis, double d_time, double rk_weight) :
    derivative{basis},
    dt{d_time},
    rkw{rk_weight}
  {}

  virtual void operator()(Sequence<Deformed_element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    // no point allocating this on the stack since face normals will usually be on the heap anyway
    Eigen::MatrixXd cartesian_normals(n_dim*n_dim, n_qpoint/row_size);
    double* cartesian_normal_ptr [n_dim];
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      auto view = cartesian_normals(Eigen::seqN(i_dim*n_dim, n_dim), Eigen::all);
      cartesian_normal_ptr[i_dim] = view.data();
      view = Eigen::VectorXd::Unit(n_dim, i_dim)*Eigen::MatrixXd::Ones(1, n_qpoint/row_size);
    }

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      auto& elem = elements[i_elem];
      double* state = elem.stage(0);
      double* rk_reference = state + n_var*n_qpoint;
      double* normals = elem.reference_level_normals();
      double* det = elem.jacobian_determinant();
      double time_rate [n_var][n_qpoint] {};
      double flux [n_var][n_dim][n_qpoint] {};
      double* face = elem.face();
      double* face_nrml [2*n_dim];
      for (int i_face = 0; i_face < 2*n_dim; ++i_face) {
        face_nrml[i_face] = elem.face_normal(i_face);
        if (!face_nrml[i_face]) {
          printf("element %p face %i\n", (void*)&elem, i_face);
          face_nrml[i_face] = cartesian_normal_ptr[i_face/2];
        }
      }
      double* tss = elem.time_step_scale();
      double* av_coef = elem.art_visc_coef();
      double nom_sz = elem.nominal_size();
      const double d_t_by_d_pos = dt/nom_sz;

      // compute gradient (with sign inverted)
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // fetch this row of data
            Eigen::Matrix<double, row_size, n_var> row_r;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_row = 0; i_row < row_size; ++i_row) {
                row_r(i_row, i_var) = state[i_var*n_qpoint + i_outer*stride*row_size + i_inner + i_row*stride];
              }
            }
            // fetch numerical face state
            Eigen::Matrix<double, 2, n_var> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                boundary_values(is_positive, i_var) = face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint];
              }
            }
            // compute gradient
            for (int j_dim = 0; j_dim < n_dim; ++j_dim)
            {
              Eigen::Matrix<double, row_size, n_var> temp_row = row_r;
              Eigen::Matrix<double, 2, n_var> temp_face = boundary_values;
              for (int i_row = 0; i_row < row_size; ++i_row) {
                temp_row(i_row, Eigen::all) *= normals[(i_dim*n_dim + j_dim)*n_qpoint + i_outer*stride*row_size + i_inner + i_row*stride];
              }
              for (int is_positive : {0, 1}) {
                temp_face(is_positive, Eigen::all) *= face_nrml[2*i_dim + is_positive][j_dim*n_qpoint/row_size + i_face_qpoint];
              }
              Eigen::Matrix<double, row_size, n_var> row_f = -derivative(temp_row, temp_face)/nom_sz;
              // add dimensional component to gradient
              for (int i_var = 0; i_var < n_var; ++i_var) {
                for (int i_row = 0; i_row < row_size; ++i_row) {
                  int offset = i_outer*stride*row_size + i_inner + i_row*stride;
                  flux[i_var][j_dim][offset] += row_f(i_row, i_var);
                }
              }
            }
            ++i_face_qpoint;
          }
        }
      }

      // compute viscous flux from gradient
      for (int i_var = 0; i_var < n_var; ++i_var) {
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          double temp_flux [n_dim] {};
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
              temp_flux[i_dim] += normals[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint]*flux[i_var][j_dim][i_qpoint];
            }
          }
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            flux[i_var][i_dim][i_qpoint] = temp_flux[i_dim]*av_coef[i_qpoint]/det[i_qpoint];
          }
        }
      }

      // compute interior term of update and write flux to faces
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // fetch this row of data
            Eigen::Matrix<double, row_size, n_var> row_f;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_row = 0; i_row < row_size; ++i_row) {
                row_f(i_row, i_var) = flux[i_var][i_dim][i_outer*stride*row_size + i_inner + i_row*stride];
              }
            }
            // differentiate flux and write to faces
            Eigen::Matrix<double, 2, n_var> boundary_values;
            Eigen::Matrix<double, row_size, n_var> row_w = -derivative.interior_term(row_f, boundary_values);
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint] = boundary_values(is_positive, i_var);
              }
            }
            // add dimensional component to update
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int i_row = 0; i_row < row_size; ++i_row) {
                int offset = i_outer*stride*row_size + i_inner + i_row*stride;
                time_rate[i_var][offset] += row_w(i_row, i_var);
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
          double updated = time_rate[i_var][i_qpoint]*d_t_by_d_pos*tss[i_qpoint]/det[i_qpoint] + state[i_dof];
          state[i_dof] = rkw*updated + (1. - rkw)*rk_reference[i_dof];
        }
      }
    }
  }
};

template<>
class Kernel_traits<Local_av0_deformed>
{
  public:
  using base_t = Kernel<Deformed_element&>;
};

}
#endif
