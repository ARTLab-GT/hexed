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

    // apply to a single element
    void operator()(const double* read, std::array<double*, 6> faces)
    {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
          auto row_r = Row_rw<Pde::n_var, row_size>::read_row(read, ind);
          Mat<2, Pde::n_var> bound = boundary*row_r;
          Row_rw<Pde::n_var, row_size>::write_bound(bound, faces, ind);
        }
      }
    }

    // apply to a sequence of elements
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
    const Pde eq;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Mat<2, row_size> boundary;
    Write_face<n_dim, row_size> write_face;
    const int stage;
    // weights for different parameters when assembling the updated state
    const double update_conv; // how much of the convective time derivative to include
    const double update_diff; // how much of the diffusive time derivative
    const double ref; // how much of the reference state to include
    const double curr; // how much of the current state to include
    double mcc;
    double mcd;

    public:
    template <typename... pde_args>
    Local(const Basis& basis, double dt, bool which_stage, pde_args... args) :
      eq(args...),
      derivative{basis},
      boundary{basis.boundary()},
      write_face{basis},
      stage{which_stage},
      update_conv{stage ? dt*basis.cancellation_convective()/basis.max_cfl_convective() : dt},
      update_diff{stage ? dt*basis.cancellation_diffusive()/basis.max_cfl_diffusive() : dt},
      ref{stage ? update_conv/dt : 0},
      curr{1 - ref},
      mcc{basis.max_cfl_convective()},
      mcd{basis.max_cfl_diffusive()}
    {
      HEXED_ASSERT(Pde::has_convection || !which_stage, "for pure diffusion use alternating time steps");
    }

    virtual void operator()(Sequence<element_t&>& elements)
    {
      // dummy face normal vector to use when deformed elements participate in cartesian connections
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
        double* time_step_scale = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [2][Pde::n_update][n_qpoint] {}; // first part contains convective time derivative, second part diffusive
        double* av_coef = elem.art_visc_coef();
        // only need the next 2 for deformed elements
        double* nrml = nullptr; // reference level normals
        double* elem_det = nullptr; // jacobian determinant
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        std::array<double*, 6> face_nrml;
        double visc_storage [n_dim][Pde::n_var][n_qpoint] {}; // for viscous PDEs, will be used to store gradients and later viscous flux
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
          // compute gradient (times jacobian determinant, cause that's easier)
          constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
          std::array<double*, 6> visc_faces;
          for (int i_face = 0; i_face < 2*n_dim; ++i_face) visc_faces[i_face] = elem.faces[i_face] + 2*(n_dim + 2)*n_qpoint/row_size;
          for (int i_dim = 0; i_dim < n_dim; ++i_dim)
          {
            for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind)
            {
              // fetch state data
              auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
              auto face_state = Row_rw<Pde::n_var, row_size>::read_bound(visc_faces, ind);
              if constexpr (element_t::is_deformed) {
                // fetch normal data
                auto row_n = Row_rw<n_dim, row_size>::read_row(nrml + i_dim*n_dim*n_qpoint, ind);
                auto face_n = Row_rw<n_dim, row_size>::read_bound(face_nrml, ind);
                // differentiate and write to temporary storage
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
                    Mat<row_size, 1> row = row_n(Eigen::all, j_dim).cwiseProduct(row_r(Eigen::all, i_var));
                    Mat<2, 1> bound = face_n(Eigen::all, j_dim).cwiseProduct(face_state(Eigen::all, i_var));
                    Row_rw<1, row_size>::write_row(derivative(row, bound)/d_pos, visc_storage[j_dim][i_var], ind, 1.);
                  }
                }
              } else {
                // differentiate and write to temporary storage
                Row_rw<Pde::n_var, row_size>::write_row(derivative(row_r, face_state)/d_pos, visc_storage[i_dim][0], ind, 0);
              }
            }
          }
          // compute viscous flux from gradient
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
          {
            // fetch normals matrix, gradient, and state
            Mat<Pde::n_var> qpoint_state;
            Mat<n_dim, Pde::n_var> qpoint_grad;
            Mat<n_dim, n_dim> qpoint_nrmls;
            qpoint_nrmls.setIdentity();
            for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
              qpoint_state(i_var) = state[i_var*n_qpoint + i_qpoint];
              for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
                qpoint_grad(i_dim, i_var) = visc_storage[i_dim][i_var][i_qpoint];
              }
            }
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                if constexpr (element_t::is_deformed) {
                  qpoint_nrmls(i_dim, j_dim) = nrml[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
                }
              }
            }
            if constexpr (element_t::is_deformed) qpoint_grad /= elem_det[i_qpoint]; // divide by the determinant to get the actual gradient (see above)
            // compute flux and write to temporary storage
            Mat<n_dim, Pde::n_var> flux = qpoint_nrmls*eq.flux_visc(qpoint_state, qpoint_grad, av_coef[i_qpoint]);
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
            Mat<row_size, Pde::n_update> flux;
            Mat<2, Pde::n_update> face_f;
            // compute convective update
            if constexpr (Pde::has_convection) {
              // fetch row data
              auto row_r = Row_rw<Pde::n_var, row_size>::read_row(state, ind);
              Mat<row_size, n_dim> row_n = Mat<row_size, 1>::Ones()*Mat<1, n_dim>::Unit(i_dim);
              if constexpr (element_t::is_deformed) {
                row_n = Row_rw<n_dim, row_size>::read_row(nrml + i_dim*n_dim*n_qpoint, ind);
              }
              // compute flux
              for (int i_row = 0; i_row < row_size; ++i_row) {
                flux(i_row, Eigen::all) = eq.flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all));
              }
              // fetch face data
              face_f = Row_rw<Pde::n_update, row_size>::read_bound(faces, ind);
              // differentiate and write to temporary storage
              Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux, face_f), time_rate[0][0], ind, 1.);
            }
            // compute viscous update
            if constexpr (Pde::is_viscous) {
              flux = Row_rw<Pde::n_update, row_size>::read_row(visc_storage[i_dim][Pde::curr_start], ind);
              face_f = boundary*flux;
              // write viscous flux to faces to enable calculation of the numerical flux
              for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
                for (int is_positive : {0, 1}) {
                  faces[ind.i_dim*2 + is_positive][(2*(n_dim + 2) + i_var)*ind.n_fqpoint + ind.i_face_qpoint()] = face_f(is_positive, i_var);
                }
              }
              Row_rw<Pde::n_update, row_size>::write_row(-derivative(flux), time_rate[1][0], ind, 1.);
            }
          }
        }

        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          Mat<Pde::n_var> s;
          for (int i_var = 0; i_var < Pde::n_var; ++i_var) s(i_var) = state[i_var*n_qpoint + i_qpoint];
          double scale = 0;
          double tss = time_step_scale[i_qpoint];
          if constexpr (Pde::has_convection) scale += eq.char_speed(s)/mcc/tss;
          if constexpr (Pde::is_viscous) scale += eq.diffusivity(s, av_coef[i_qpoint])/mcd/tss/tss;
          scale *= n_dim;
          for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
            if constexpr (Pde::has_convection) time_rate[0][i_var][i_qpoint] /= scale;
            if constexpr (Pde::is_viscous)     time_rate[1][i_var][i_qpoint] /= scale;
          }
        }

        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var)
        {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* ref_state = state + (Pde::ref_start + i_var)*n_qpoint;
          double* visc_state = nullptr; // only need this for viscous, of course
          if constexpr (Pde::is_viscous) visc_state = state + (Pde::visc_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            if constexpr (Pde::has_convection) {
              double u = update_conv*time_rate[0][i_var][i_qpoint];
              if constexpr (Pde::is_viscous) {
                u += update_diff*time_rate[1][i_var][i_qpoint];
                // for stage 1, the updates for convection and diffusion are different,
                // so we have to mix in some of the diffisive update from stage 0
                // to achieve proper cancellation
                if (stage) u += (update_conv - update_diff)*visc_state[i_qpoint];
                else visc_state[i_qpoint] = time_rate[1][i_var][i_qpoint]; // for stage 0, record the diffusive update to allow the aforementioned
              }
              u /= det*d_pos;
              curr_state[i_qpoint] = u + curr*curr_state[i_qpoint] + ref*ref_state[i_qpoint];
            } else {
              curr_state[i_qpoint] += update_diff*time_rate[1][i_var][i_qpoint]/det/d_pos;
            }
          }
        }
        // write updated state to face storage.
        // For viscous, don't bother since we still have to add the numerical flux term
        if constexpr (!Pde::is_viscous) write_face(state, elem.faces);
      }
    }
  };

  // account for the difference between the real and numerical viscous fluxes to enforce conservation
  template <int n_dim, int row_size>
  class Reconcile_ldg_flux : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    Pde eq;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    int stage;
    double update;
    double mcc;
    double mcd;

    public:
    template <typename... pde_args>
    Reconcile_ldg_flux(const Basis& basis, double dt, int which_stage, pde_args... args) :
      eq(args...),
      derivative{basis},
      write_face{basis},
      stage{which_stage},
      update{stage ? dt*basis.cancellation_diffusive()/basis.max_cfl_diffusive() : dt},
      mcc{basis.max_cfl_convective()},
      mcd{basis.max_cfl_diffusive()}
    {
      HEXED_ASSERT(Pde::has_convection || !which_stage, "for pure diffusion use alternating time steps");
    }

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* state = elem.stage(0);
        std::array<double*, 6> faces = elem.faces;
        for (double*& face : faces) face += (2*(n_dim + 2) + Pde::curr_start)*n_qpoint/row_size;
        double* time_step_scale = elem.time_step_scale();
        double* av_coef = elem.art_visc_coef();
        double d_pos = elem.nominal_size();
        double time_rate [Pde::n_update][n_qpoint] {};
        double* elem_det = nullptr;
        if constexpr (element_t::is_deformed) {
          elem_det = elem.jacobian_determinant();
        }
        // fetch the flux difference and lift it to the interior temporary storage
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(faces, ind);
            Row_rw<Pde::n_update, row_size>::write_row(-derivative.boundary_term(bound_f), time_rate[0], ind, 1.);
          }
        }

        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          Mat<Pde::n_var> s;
          for (int i_var = 0; i_var < Pde::n_var; ++i_var) s(i_var) = state[i_var*n_qpoint + i_qpoint];
          double scale = 0;
          double tss = time_step_scale[i_qpoint];
          if constexpr (Pde::has_convection) scale += eq.char_speed(s)/mcc/tss;
          if constexpr (Pde::is_viscous) scale += eq.diffusivity(s, av_coef[i_qpoint])/mcd/tss/tss;
          scale *= n_dim;
          for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
            time_rate[i_var][i_qpoint] /= scale;
          }
        }

        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* visc_state = state + (Pde::visc_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            curr_state[i_qpoint] += update*time_rate[i_var][i_qpoint]/d_pos/det;
            if constexpr (Pde::has_convection) if (!stage) visc_state[i_qpoint] += time_rate[i_var][i_qpoint]; // add face flux correction to viscous update to be used in stage 1
          }
        }
        // *now* we can extrapolate state to faces
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
    const Pde eq;
    static constexpr int n_fqpoint = custom_math::pow(row_size, n_dim - 1);

    public:
    Neighbor() : eq{} {}

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
          auto flux = eq.flux_num(state, nrml);
          // write flux to temporary storage
          for (int i_side = 0; i_side < 2; ++i_side) {
            for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
              face[i_side][(i_var + Pde::curr_start)*n_fqpoint + i_qpoint] = sign[i_side]*flux(i_var);
            }
            // compute average face state for LDG scheme
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
          // write flux
          double* f = con.state() + i_side*n_fqpoint*(n_dim + 2);
          int offset = Pde::curr_start*n_fqpoint;
          for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
            f[offset + i_dof] = face[i_side][offset + i_dof];
          }
          // write average face state
          if constexpr (Pde::is_viscous) {
            for (int i_dof = 0; i_dof < Pde::n_var*n_fqpoint; ++i_dof) {
              f[2*(n_dim + 2)*n_fqpoint + i_dof] = face[2 + i_side][i_dof];
            }
          }
        }
      }
    }
  };

  // compute the difference between the numerical (average) viscous flux and
  // the viscous flux on each face in preparation for reconciling the different face fluxes
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
      }
    }
  };

  // compute the maximum stable time step
  template <int n_dim, int row_size>
  class Max_dt : public Kernel<element_t&, double>
  {
    using Pde = Pde_templ<n_dim>;
    const Pde eq;
    double max_cfl_c;
    double max_cfl_d;

    public:
    Max_dt(const Basis& basis) :
      eq{},
      max_cfl_c{basis.max_cfl_convective()},
      max_cfl_d{basis.max_cfl_diffusive()}
    {}

    virtual double operator()(Sequence<element_t&>& elements)
    {
      constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
      // compute the maximum stable time step for all elements and take the minimum
      double dt = std::numeric_limits<double>::max();
      #pragma omp parallel for reduction(min:dt)
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        element_t& elem {elements[i_elem]};
        double* state = elem.stage(0);
        double* art_visc = elem.art_visc_coef();
        double* tss = elem.time_step_scale();
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        double* ref_nrml = nullptr;
        double* jac_det = nullptr;
        #pragma GCC diagnostic pop
        if constexpr (element_t::is_deformed) {
          ref_nrml = elem.reference_level_normals();
          jac_det = elem.jacobian_determinant();
        }
        // minimum spatial scales (determined by Jacobian) for convection and diffusion
        double min_scale_conv = std::numeric_limits<double>::max(); // O(h)
        double min_scale_diff = std::numeric_limits<double>::max(); // O(h^2)
        // effective speeds of information propagation
        double max_qpoint_speed = 0.; // maximum characteristic speed
        double max_qpoint_diff = 0.; // maximum effective diffusivity
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
        {
          // fetch state variables
          Mat<Pde::n_var> qpoint_state;
          for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
            qpoint_state(i_var) = state[i_var*n_qpoint + i_qpoint];
          }
          // compute speeds
          if constexpr (Pde::has_convection) max_qpoint_speed = std::max(max_qpoint_speed, eq.char_speed(qpoint_state));
          if constexpr (Pde::is_viscous) max_qpoint_diff = std::max(max_qpoint_diff, eq.diffusivity(qpoint_state, art_visc[i_qpoint]));
          // get time step scale (which for pure diffusion equations will generally be squared
          double local_tss = custom_math::pow(tss[i_qpoint], Pde::tss_pow);
          // compute spatial scales
          double local_size = elem.nominal_size();
          if constexpr (element_t::is_deformed) {
            double norm_sum = 0.;
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              double norm_sq = 0.;
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                double coef = ref_nrml[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
                norm_sq += coef*coef;
              }
              norm_sum += std::sqrt(norm_sq);
            }
            local_size *= n_dim*jac_det[i_qpoint]/norm_sum; // for deformed elements this is a essentially a measure of the amount of stretching in each dimension
          }
          // take minima
          min_scale_conv = std::min(min_scale_conv, local_size/local_tss);
          min_scale_diff = std::min(min_scale_diff, local_size*local_size/local_tss);
        }
        // combine convective and diffusive time steps in a way that is less than either,
        // but asymptotes to the dominant one (if any)
        dt = std::min(dt, 1./n_dim/(max_qpoint_speed/max_cfl_c/min_scale_conv + max_qpoint_diff/max_cfl_d/min_scale_diff));
      }
      return dt;
    }
  };
};

}
#endif
