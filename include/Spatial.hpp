#ifndef HEXED_SPATIAL_HPP_
#define HEXED_SPATIAL_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Write_face.hpp"
#include "Row_rw.hpp"

namespace hexed
{

template <typename element_t, template<int> typename Pde_templ>
class Spatial
{
  public:
  Spatial() = delete;

  template <int n_dim, int row_size>
  class Local : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    using Rrw = Row_rw<Pde::n_var, row_size>;
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    double update;
    double curr;
    double ref;
    const double heat_rat;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);

    public:
    Local(const Basis& basis,
          double update_coef, double current_coef, double reference_coef,
          double heat_ratio=1.4) :
      derivative{basis},
      write_face{basis},
      update{update_coef},
      curr{current_coef},
      ref{reference_coef},
      heat_rat{heat_ratio}
    {}

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* active_state = elem.stage(0);
        double* reference_state = active_state + Pde::n_var*n_qpoint;
        double* face = elem.face();
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [Pde::n_var][n_qpoint] {};

        // compute update
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            auto row_r = Rrw::read_row(active_state, ind);
            Mat<row_size, n_dim> row_n = Mat<row_size, 1>::Ones()*Mat<1, n_dim>::Unit(i_dim);
            if constexpr (element_t::is_deformed) {
              row_n = Row_rw<n_dim, row_size>::read_row(elem.reference_level_normals() + i_dim*n_dim*n_qpoint, ind);
            }
            Mat<row_size, Pde::n_var> flux;
            for (int i_row = 0; i_row < row_size; ++i_row) {
              flux(i_row, Eigen::all) = Pde::flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all), i_dim);
            }
            Rrw::write_row(-derivative(flux, Rrw::read_bound(face, ind)), time_rate[0], ind, 1.);
          }
        }

        // write the updated solution
        double* elem_det = nullptr;
        if constexpr (element_t::is_deformed) elem_det = elem.jacobian_determinant();
        for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            const int i_dof = i_var*n_qpoint + i_qpoint;
            active_state[i_dof] = update*time_rate[i_var][i_qpoint]/d_pos*tss[i_qpoint]/det
                                  + curr*active_state[i_dof]
                                  + ref*reference_state[i_dof];
          }
        }
        write_face(active_state, face);
      }
    }
  };
};

}
#endif
