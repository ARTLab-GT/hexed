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

// class to contain all the kernels in a scope
// parameterized by the type of element (deformed/cartesian) and the PDE
template <typename element_t, template<int> typename Pde_templ>
class Spatial
{
  public:
  // this class has no non-static data members, so there is no reason to construct it
  Spatial() = delete;
  Spatial(const Spatial&) = delete;
  Spatial(Spatial&&) = delete;

  // extrapolates the values in the interior of an element to the faces
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

    void operator()(const double* read, std::array<double*, 6> faces)
    {
      //read += Pde::curr_start*custom_math::pow(row_size, n_dim);
      //face += Pde::curr_start*custom_math::pow(row_size, n_dim - 1);
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
          auto row_r = Row_rw<Pde::n_var, row_size>::read_row(read, ind);
          Mat<2, Pde::n_var> bound = boundary*row_r;
          Row_rw<Pde::n_var, row_size>::write_bound(bound, faces, ind);
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
        auto faces = elem.faces;
        operator()(read, face);
        operator()(read, faces);
      }
    }
  };

  // performs the update to the element state after the shared numerical flux
  // has been computed.
  // In other words, performs the part of the algorithm that only uses local data
  // (or in other words does not directly depend on neighboring elements).
  // Note that face data is updated to reflect the updated interior state.
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

    // note: `.*_coef` are the coefficients for multistage time integration.
    // The final result is
    // `update_coef`*[residual] + `current_coef`*[current state] + `reference_coef`*[reference state]
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

        // compute residual
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            // fetch row data
            auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
            Mat<row_size, n_dim> row_n = Mat<row_size, 1>::Ones()*Mat<1, n_dim>::Unit(i_dim);
            if constexpr (element_t::is_deformed) {
              row_n = Row_rw<n_dim, row_size>::read_row(elem.reference_level_normals() + i_dim*n_dim*n_qpoint, ind);
            }
            // compute flux
            Mat<row_size, Pde::n_update> flux;
            for (int i_row = 0; i_row < row_size; ++i_row) {
              flux(i_row, Eigen::all) = Pde::flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all));
            }
            // fetch boundary data
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(face + Pde::curr_start*n_qpoint/row_size, ind);
            // differentiate and write to temporary storage
            Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux, bound_f), time_rate[0], ind, 1.);
          }
        }

        // write update to interior
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
        // write updated state to face storage
        write_face(state, face);
        write_face(state, elem.faces);
      }
    }
  };

  // Computes the shared numerical flux at the element interfaces.
  // Requires that the state has been written to the face storage
  // and replaces the state of both faces with the computed flux.
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
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable" // otherwise `'face_nrml' set but not used` (not sure why)
        auto& con = connections[i_con];
        auto dir = con.direction();
        double face [2][(n_dim + 2)*n_fqpoint]; // copying face data to temporary stack storage improves efficiency
        double face_nrml [n_dim*n_fqpoint]; // only set for deformed
        int sign [2] {1, 1}; // records whether the normal vector on each side needs to be flipped to obey sign convention
        // fetch face data
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.face(i_side);
          for (int i_dof = 0; i_dof < Pde::n_var*n_fqpoint; ++i_dof) {
            face[i_side][i_dof] = f[i_dof];
          }
        }
        Mat<n_dim> nrml; // will be used in loop to contain reference level normal
        Face_permutation<n_dim, row_size> perm(dir, face[1]); // only used for deformed
        if constexpr (element_t::is_deformed) {
          perm.match_faces(); // if order of quadrature points on both faces does not match, reorder face 1 to match face 0
          for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
          double* n = con.normal();
          for (int i_dof = 0; i_dof < n_dim*n_fqpoint; ++i_dof) {
            face_nrml[i_dof] = n[i_dof];
          }
        } else {
          nrml.setUnit(dir.i_dim); // normal vector is trivial for cartesian
        }
        // compute flux
        for (int i_qpoint = 0; i_qpoint < n_fqpoint; ++i_qpoint)
        {
          // fetch data
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
          // compute flux
          auto flux = Pde::flux_num(state, nrml);
          // write flux to temporary storage
          for (int i_side = 0; i_side < 2; ++i_side) {
            for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
              face[i_side][(i_var + Pde::curr_start)*n_fqpoint + i_qpoint] = sign[i_side]*flux(i_var);
            }
          }
        }
        if constexpr (element_t::is_deformed) perm.restore(); // restore data of face 1 to original order
        // write data to actual face storage on heap
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
