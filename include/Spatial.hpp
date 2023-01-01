#ifndef HEXED_SPATIAL_HPP_
#define HEXED_SPATIAL_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Row_rw.hpp"

namespace hexed
{

template <typename element_t, template<int> typename Pde_templ>
class Spatial
{
  public:
  Spatial() = delete;
  Spatial(const Spatial&) = delete;
  Spatial(Spatial&&) = delete;

  template <int n_dim, int row_size>
  class Write_face : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    const Eigen::Matrix<double, 2, row_size> boundary;

    public:
    Write_face(const Basis& basis) : boundary{basis.boundary()} {}

    void operator()(const double* read, double* face)
    {
      //read += Pde::curr_start*custom_math::pow(row_size, n_dim);
      //face += Pde::curr_start*custom_math::pow(row_size, n_dim - 1);
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
          auto row_r = Row_rw<Pde::n_var, row_size>::read_row(read, ind);
          Mat<2, Pde::n_var> bound = boundary*row_r;
          Row_rw<Pde::n_var, row_size>::write_bound(bound, face, ind);
        }
      }
    }

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
        element_t& elem {elements[i_elem]};
        double* read = elem.stage(0);
        double* face = elem.face();
        operator()(read, face);
      }
    }
  };

  template <int n_dim, int row_size>
  class Local : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    double update;
    double curr;
    double ref;
    const double heat_rat;

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
        double* state = elem.stage(0);
        double* face = elem.face();
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [Pde::n_update][n_qpoint] {};

        // compute update
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
            Mat<row_size, n_dim> row_n = Mat<row_size, 1>::Ones()*Mat<1, n_dim>::Unit(i_dim);
            if constexpr (element_t::is_deformed) {
              row_n = Row_rw<n_dim, row_size>::read_row(elem.reference_level_normals() + i_dim*n_dim*n_qpoint, ind);
            }
            Mat<row_size, Pde::n_update> flux;
            for (int i_row = 0; i_row < row_size; ++i_row) {
              flux(i_row, Eigen::all) = Pde::flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all));
            }
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(face + Pde::curr_start*n_qpoint/row_size, ind);
            Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux, bound_f), time_rate[0], ind, 1.);
          }
        }

        // write the updated solution
        double* elem_det = nullptr;
        if constexpr (element_t::is_deformed) elem_det = elem.jacobian_determinant();
        for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* ref_state = state + (Pde::ref_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            curr_state[i_qpoint] = update*time_rate[i_var][i_qpoint]/d_pos*tss[i_qpoint]/det
                                   + curr*curr_state[i_qpoint] + ref*ref_state[i_qpoint];
          }
        }
        write_face(state, face);
      }
    }
  };

  template <int n_dim, int row_size>
  class Neighbor : public Kernel<Face_connection<element_t>&>
  {
    using Pde = Pde_templ<n_dim>;
    static constexpr int n_fqpoint = custom_math::pow(row_size, n_dim - 1);
    public:
    virtual void operator()(Sequence<Face_connection<element_t>&>& connections)
    {
      #pragma omp parallel for
      for (int i_con = 0; i_con < connections.size(); ++i_con)
      {
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable" // otherwise `'face_nrml' set but not used`
        auto& con = connections[i_con];
        auto dir = con.direction();
        double face [2][Pde::n_var*n_fqpoint];
        double face_nrml [n_dim*n_fqpoint];
        int sign [2] {1, 1};
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.face(i_side);
          for (int i_dof = 0; i_dof < Pde::n_var*n_fqpoint; ++i_dof) {
            face[i_side][i_dof] = f[i_dof];
          }
        }
        Mat<n_dim> nrml;
        Face_permutation<n_dim, row_size> perm(dir, face[1]);
        if constexpr (element_t::is_deformed) {
          perm.match_faces();
          for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
          double* n = con.normal();
          for (int i_dof = 0; i_dof < n_dim*n_fqpoint; ++i_dof) {
            face_nrml[i_dof] = n[i_dof];
          }
        } else {
          nrml.setUnit(dir.i_dim);
        }
        for (int i_qpoint = 0; i_qpoint < n_fqpoint; ++i_qpoint)
        {
          if constexpr (element_t::is_deformed) {
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              nrml(i_dim) = sign[0]*face_nrml[i_dim*n_fqpoint + i_qpoint];
            }
          }
          Mat<Pde::n_var, 2> state;
          for (int i_side = 0; i_side < 2; ++i_side) {
            for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
              state(i_var, i_side) = face[i_side][i_var*n_fqpoint + i_qpoint];
            }
          }
          auto flux = Pde::flux_num(state, nrml);
          for (int i_side = 0; i_side < 2; ++i_side) {
            for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
              face[i_side][(i_var + Pde::curr_start)*n_fqpoint + i_qpoint] = sign[i_side]*flux(i_var);
            }
          }
        }
        if constexpr (element_t::is_deformed) perm.restore();
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.face(i_side);
          int offset = Pde::curr_start*n_fqpoint;
          for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
            f[offset + i_dof] = face[i_side][offset + i_dof];
          }
        }
        #pragma GCC diagnostic pop
      }
    }
  };
};

}
#endif
