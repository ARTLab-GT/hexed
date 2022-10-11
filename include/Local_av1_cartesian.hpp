#ifndef HEXED_LOCAL_AV1_CARTESIAN_HPP_
#define HEXED_LOCAL_AV1_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Write_face.hpp"

namespace hexed
{

/*
 * Computes the local update for the first phase of the artificial viscosity operator (for Cartesian elements).
 * Requires that faces contain the shared LDG state.
 */
template <int n_dim, int row_size>
class Local_av1_cartesian : public Kernel<Element&>
{
  Derivative<row_size> derivative;
  Write_face<n_dim, row_size> write_face;
  double dt;
  bool use_coef;

  public:
  Local_av1_cartesian(const Basis& basis, double d_time, bool use_av_coef = true) :
    derivative{basis},
    write_face{basis},
    dt{d_time},
    use_coef{use_av_coef}
  {}

  virtual void operator()(Sequence<Element&>& elements)
  {
    constexpr int n_var = n_dim + 2;
    constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    constexpr int n_face_dof = n_var*n_qpoint/row_size;

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
    {
      auto& elem = elements[i_elem];
      double* state = elem.stage(0);
      double time_rate [n_var][n_qpoint] {};
      double* face = elem.face();
      double* tss = elem.time_step_scale();
      double nom_sz = elem.nominal_size();
      const double d_t_by_d_pos = dt/nom_sz;

      // Compute update
      for (int stride = n_qpoint/row_size, n_rows = 1, i_dim = 0; n_rows < n_qpoint;
           stride /= row_size, n_rows *= row_size, ++i_dim)
      {
        int i_face_qpoint {0};
        for (int i_outer = 0; i_outer < n_rows; ++i_outer)
        {
          for (int i_inner = 0; i_inner < stride; ++i_inner)
          {
            // fetch numerical face flux
            Eigen::Matrix<double, 2, n_var> boundary_values;
            for (int i_var = 0; i_var < n_var; ++i_var) {
              for (int is_positive : {0, 1}) {
                boundary_values(is_positive, i_var) = face[(i_dim*2 + is_positive)*n_face_dof + i_var*n_qpoint/row_size + i_face_qpoint];
              }
            }
            Eigen::Matrix<double, row_size, n_var> row_w = -derivative.boundary_term(boundary_values);
            // Add dimensional component to update
            for (int i_var = 0; i_var < n_var; ++i_var) {
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
          state[i_dof] += time_rate[i_var][i_qpoint]*d_t_by_d_pos*tss[i_qpoint]
                          *(use_coef ? 1 : tss[i_qpoint]);
        }
      }
      write_face(state, face);
    }
  }
};

template<>
class Kernel_traits<Local_av1_cartesian>
{
  public:
  using base_t = Kernel<Element&>;
};

}
#endif
