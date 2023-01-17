#ifndef HEXED_SPATIAL_HPP_
#define HEXED_SPATIAL_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Row_rw.hpp"
#include "connection.hpp"
#include "Face_permutation.hpp"

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
        auto faces = elem.faces;
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
    Mat<2, row_size> boundary;
    Write_face<n_dim, row_size> write_face;
    const int stage;
    const double update_conv;
    const double update_diff;
    const double ref;
    const double curr;
    const double heat_rat;

    // note: `.*_coef` are the coefficients for multistage time integration.
    // The final result is
    // `update_coef`*[residual] + `current_coef`*[current state] + `reference_coef`*[reference state]
    public:
    Local(const Basis& basis, double dt, bool which_stage, double heat_ratio=1.4) :
      derivative{basis},
      boundary{basis.boundary()},
      write_face{basis},
      stage{which_stage},
      update_conv{stage ? dt*basis.cancellation_convective()/basis.max_cfl_convective() : dt},
      update_diff{stage ? dt*basis.cancellation_diffusive()/basis.max_cfl_diffusive() : dt},
      ref{stage ? update_conv/dt : 0},
      curr{1 - ref},
      heat_rat{heat_ratio}
    {}

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
      double cartesian_normal [n_dim][n_dim][n_qpoint/row_size];
      #pragma GCC diagnostic pop
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          for (int i_qpoint = 0; i_qpoint < n_qpoint/row_size; ++i_qpoint) {
            cartesian_normal[i_dim][j_dim][i_qpoint] = i_dim == j_dim;
          }
        }
      }

      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* state = elem.stage(0);
        std::array<double*, 6> faces = elem.faces;
        for (double*& face : faces) face += Pde::curr_start*n_qpoint/row_size;
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [2][Pde::n_update][n_qpoint] {};
        double* av_coef = elem.art_visc_coef();
        double* nrml = nullptr;
        double* elem_det = nullptr;
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        std::array<double*, 6> face_nrml;
        double visc_storage [n_dim][Pde::n_var][n_qpoint] {};
        #pragma GCC diagnostic pop
        if constexpr (element_t::is_deformed) {
          elem_det = elem.jacobian_determinant();
          nrml = elem.reference_level_normals();
          for (int i_face = 0; i_face < 2*n_dim; ++i_face) {
            face_nrml[i_face] = elem.face_normal(i_face);
            if (!face_nrml[i_face]) face_nrml[i_face] = cartesian_normal[i_face/2][0];
          }
        }

        if constexpr (Pde::is_viscous)
        {
          // compute gradient
          constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
          std::array<double*, 6> visc_faces;
          for (int i_face = 0; i_face < 2*n_dim; ++i_face) visc_faces[i_face] = elem.faces[i_face] + 2*(n_dim + 2)*n_qpoint/row_size;
          for (int i_dim = 0; i_dim < n_dim; ++i_dim)
          {
            for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind)
            {
              // fetch state data
              auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
              auto bound_state = Row_rw<Pde::n_var, row_size>::read_bound(visc_faces, ind);
              if constexpr (element_t::is_deformed) {
                // fetch normal data
                auto row_n = Row_rw<n_dim, row_size>::read_row(nrml + i_dim*n_dim*n_qpoint, ind);
                auto bound_nrml = Row_rw<n_dim, row_size>::read_bound(face_nrml, ind);
                // differentiate and write to temporary storage
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
                    Mat<row_size, 1> row = row_n(Eigen::all, j_dim).cwiseProduct(row_r(Eigen::all, i_var));
                    Mat<2, 1> bound = bound_nrml(Eigen::all, j_dim).cwiseProduct(bound_state(Eigen::all, i_var));
                    Row_rw<1, row_size>::write_row(derivative(row, bound)/d_pos, visc_storage[j_dim][i_var], ind, 1.);
                  }
                }
              } else {
                Row_rw<Pde::n_var, row_size>::write_row(derivative(row_r, bound_state)/d_pos, visc_storage[i_dim][0], ind, 0);
              }
            }
          }
          // compute viscous flux from gradient
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
          {
            Mat<n_dim, Pde::n_var> qpoint_grad;
            Mat<n_dim, n_dim> qpoint_nrmls;
            qpoint_nrmls.setIdentity();
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int i_var = 0; i_var < Pde::n_var; ++i_var) qpoint_grad(i_dim, i_var) = visc_storage[i_dim][i_var][i_qpoint];
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                if constexpr (element_t::is_deformed) {
                  qpoint_nrmls(i_dim, j_dim) = nrml[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
                }
              }
            }
            if constexpr (element_t::is_deformed) qpoint_grad /= elem_det[i_qpoint];
            auto flux = Pde::flux_visc(qpoint_grad, qpoint_nrmls, av_coef[i_qpoint]);
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
                visc_storage[i_dim][Pde::curr_start + i_var][i_qpoint] = flux(i_dim, i_var);
              }
            }
          }
        }

        // compute residual
        for (int i_dim = 0; i_dim < n_dim; ++i_dim)
        {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind)
          {
            // compute convective update
            // fetch row data
            auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
            Mat<row_size, n_dim> row_n = Mat<row_size, 1>::Ones()*Mat<1, n_dim>::Unit(i_dim);
            if constexpr (element_t::is_deformed) {
              row_n = Row_rw<n_dim, row_size>::read_row(nrml + i_dim*n_dim*n_qpoint, ind);
            }
            // compute flux
            Mat<row_size, Pde::n_update> flux;
            for (int i_row = 0; i_row < row_size; ++i_row) {
              flux(i_row, Eigen::all) = Pde::flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all));
            }
            // fetch boundary data
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(faces, ind);
            // differentiate and write to temporary storage
            Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux, bound_f), time_rate[0][0], ind, 1.);
            // compute viscous update
            if constexpr (Pde::is_viscous) {
              flux = Row_rw<Pde::n_update, row_size>::read_row(visc_storage[i_dim][Pde::curr_start], ind);
              bound_f = boundary*flux;
              for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
                for (int is_positive : {0, 1}) {
                  faces[ind.i_dim*2 + is_positive][(2*(n_dim + 2) + i_var)*ind.n_fqpoint + ind.i_face_qpoint()] = bound_f(is_positive, i_var);
                }
              }
              Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux), time_rate[1][0], ind, 1.);
            }
          }
        }

        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var)
        {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* ref_state = state + (Pde::ref_start + i_var)*n_qpoint;
          double* visc_state = nullptr;
          if constexpr (Pde::is_viscous) visc_state = state + (Pde::visc_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            double u = update_conv*time_rate[0][i_var][i_qpoint];
            if constexpr (Pde::is_viscous) {
              u += update_diff*time_rate[1][i_var][i_qpoint];
              if (stage) u += (update_conv - update_diff)*visc_state[i_qpoint];
              else visc_state[i_qpoint] = time_rate[1][i_var][i_qpoint];
            }
            u *= custom_math::pow(tss[i_qpoint], Pde::tss_pow)/det/d_pos;
            curr_state[i_qpoint] = u + curr*curr_state[i_qpoint] + ref*ref_state[i_qpoint];
          }
        }
        // write updated state to face storage
        if constexpr (!Pde::is_viscous) write_face(state, elem.faces);
      }
    }
  };

  template <int n_dim, int row_size>
  class Reconcile_ldg_flux : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    int stage;
    double update;

    public:
    Reconcile_ldg_flux(const Basis& basis, double dt, int which_stage) :
      derivative{basis},
      write_face{basis},
      stage{which_stage},
      update{stage ? dt*basis.cancellation_diffusive()/basis.max_cfl_diffusive() : dt}
    {}

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* state = elem.stage(0);
        std::array<double*, 6> faces = elem.faces;
        for (double*& face : faces) face += (2*(n_dim + 2) + Pde::curr_start)*n_qpoint/row_size;
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [Pde::n_update][n_qpoint] {};
        double* elem_det = nullptr;
        if constexpr (element_t::is_deformed) {
          elem_det = elem.jacobian_determinant();
        }
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(faces, ind);
            Row_rw<Pde::n_update, row_size>::write_row(-derivative.boundary_term(bound_f), time_rate[0], ind, 1.);
          }
        }
        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* visc_state = state + (Pde::visc_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            curr_state[i_qpoint] += update*time_rate[i_var][i_qpoint]/d_pos*tss[i_qpoint]/det;
            if (!stage) visc_state[i_qpoint] += time_rate[i_var][i_qpoint];
          }
        }
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
        auto& con = connections[i_con];
        auto dir = con.direction();
        double face [4][(n_dim + 2)*n_fqpoint]; // copying face data to temporary stack storage improves efficiency
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        double face_nrml [n_dim*n_fqpoint]; // only set for deformed
        #pragma GCC diagnostic pop
        int sign [2] {1, 1}; // records whether the normal vector on each side needs to be flipped to obey sign convention
        // fetch face data
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state() + i_side*n_fqpoint*(n_dim + 2);
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
            if constexpr (Pde::is_viscous) {
              for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
                face[2 + i_side][i_var*n_fqpoint + i_qpoint] = .5*(state(i_var, 0) + state(i_var, 1));
              }
            }
          }
        }
        if constexpr (element_t::is_deformed) {
          perm.restore(); // restore data of face 1 to original order
          if constexpr (Pde::is_viscous) Face_permutation<n_dim, row_size>(dir, face[3]).restore();
        }
        // write data to actual face storage on heap
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state() + i_side*n_fqpoint*(n_dim + 2);
          int offset = Pde::curr_start*n_fqpoint;
          for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
            f[offset + i_dof] = face[i_side][offset + i_dof];
          }
          if constexpr (Pde::is_viscous) {
            for (int i_dof = 0; i_dof < Pde::n_var*n_fqpoint; ++i_dof) {
              f[2*(n_dim + 2)*n_fqpoint + i_dof] = face[2 + i_side][i_dof];
            }
          }
        }
        #pragma GCC diagnostic pop
      }
    }
  };

  template <int n_dim, int row_size>
  class Neighbor_reconcile : public Kernel<Face_connection<element_t>&>
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
        int sign [2] {1, 1}; // records whether the normal vector on each side needs to be flipped to obey sign convention
        // fetch face data
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state() + (2 + i_side)*n_fqpoint*(n_dim + 2);
          for (int i_dof = 0; i_dof < Pde::n_var*n_fqpoint; ++i_dof) {
            face[i_side][i_dof] = f[i_dof];
          }
        }
        Face_permutation<n_dim, row_size> perm(dir, face[1]); // only used for deformed
        if constexpr (element_t::is_deformed) {
          perm.match_faces(); // if order of quadrature points on both faces does not match, reorder face 1 to match face 0
          for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
        }
        // compute flux
        for (int i_qpoint = 0; i_qpoint < n_fqpoint; ++i_qpoint) {
          // fetch data
          for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
            double avg = 0;
            for (int i_side = 0; i_side < 2; ++i_side) {
              avg += .5*sign[i_side]*face[i_side][i_var*n_fqpoint + i_qpoint];
            }
            for (int i_side = 0; i_side < 2; ++i_side) {
              double& f = face[i_side][i_var*n_fqpoint + i_qpoint];
              f = sign[i_side]*avg - f;
            }
          }
        }
        if constexpr (element_t::is_deformed) perm.restore(); // restore data of face 1 to original order
        // write data to actual face storage on heap
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state() + (2 + i_side)*n_fqpoint*(n_dim + 2);
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
