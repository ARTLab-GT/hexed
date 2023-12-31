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

//! class to contain all the kernels in a scope
//! parameterized by the type of element (deformed/cartesian) and the PDE
template <typename element_t, template<int> typename Pde_templ>
class Spatial
{
  public:
  //! this class has no non-static data members, so there is no reason to construct it
  Spatial() = delete;
  Spatial(const Spatial&) = delete;
  Spatial(Spatial&&) = delete;

  //! extrapolates the values in the interior of an element to the faces
  template <int n_dim, int row_size>
  class Write_face : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    const Eigen::Matrix<double, 2, row_size> boundary;

    public:
    Write_face(const Basis& basis) : boundary{basis.boundary()} {}

    //! apply to a single element
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

    //! apply to a sequence of elements
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

  /*! \brief Performs the update to the element state after the shared numerical flux has been computed.
   * \details In other words, performs the part of the algorithm that only uses local data
   * (or in other words does not directly depend on neighboring elements).
   * Note that face data is updated to reflect the updated interior state.
   */
  template <int n_dim, int row_size>
  class Local : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    const Pde eq;
    static constexpr int n_qpoint = math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Mat<2, row_size> boundary;
    Write_face<n_dim, row_size> write_face;
    Mat<row_size, row_size> filter;
    const double _update;
    const int _stage;
    const bool _compute_residual;
    const bool _use_filter;
    // weights for different parameters when assembling the updated state

    public:
    template <typename... pde_args>
    Local(const Basis& basis, double dt, bool stage, bool compute_residual, bool use_filter, pde_args... args) :
      eq(args...),
      derivative{basis},
      boundary{basis.boundary()},
      write_face{basis},
      filter{basis.filter()},
      _update{stage ? dt*basis.step_ratio() : dt},
      _stage{stage},
      _compute_residual{compute_residual},
      _use_filter{use_filter}
    {
      HEXED_ASSERT(Pde::has_convection || !_stage, "two-stage stabilization is not applicable to diffusion equations");
      HEXED_ASSERT(!(_stage && _compute_residual), "residual calculation is a single-stage operation");
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
        double* tss = elem.time_step_scale();
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

        if constexpr (Pde::is_viscous) if (!_stage)
        {
          // compute gradient (times jacobian determinant, cause that's easier)
          constexpr int n_qpoint = math::pow(row_size, n_dim);
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
            if constexpr (Pde::is_viscous) if (!_stage) {
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

        // apply mode filtering
        if (_use_filter) {
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
              Mat<row_size, 2*Pde::n_update> row_r = Row_rw<2*Pde::n_update, row_size>::read_row(time_rate[0][0], ind);
              row_r = filter*row_r;
              Row_rw<2*Pde::n_update, row_size>::write_row(row_r, time_rate[0][0], ind, 0.);
            }
          }
        }

        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var)
        {
          double* curr_state = state + (Pde::curr_start + i_var)*n_qpoint;
          double* ref_state = state + (Pde::ref_start + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            double u = time_rate[0][i_var][i_qpoint];
            if (_stage) {
              u -= ref_state[i_qpoint];
            } else {
              if constexpr (Pde::has_convection) ref_state[i_qpoint] = u;
              if constexpr (Pde::is_viscous) u += time_rate[1][i_var][i_qpoint];
            }
            u *= _update*tss[i_qpoint]/det/d_pos;
            if (_compute_residual) ref_state[i_qpoint] = u;
            else curr_state[i_qpoint] += u;
          }
        }
        // write updated state to face storage.
        // For viscous, don't bother since we still have to add the numerical flux term
        if (!Pde::is_viscous || _stage) write_face(state, elem.faces);
      }
    }
  };

  //! account for the difference between the real and numerical viscous fluxes to enforce conservation
  template <int n_dim, int row_size>
  class Reconcile_ldg_flux : public Kernel<element_t&>
  {
    using Pde = Pde_templ<n_dim>;
    static constexpr int n_qpoint = math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    Mat<row_size, row_size> filter;
    double _update;
    int _stage;
    bool _compute_residual;
    bool _use_filter;

    public:
    Reconcile_ldg_flux(const Basis& basis, double dt, int which_stage, bool compute_residual, bool use_filter) :
      derivative{basis},
      write_face{basis},
      filter{basis.filter()},
      _update{dt},
      _stage{which_stage},
      _compute_residual{compute_residual},
      _use_filter{use_filter}
    {
      HEXED_ASSERT(Pde::has_convection || !_stage, "two-stage stabilization is not applicable to diffusion equations");
      HEXED_ASSERT(Pde::has_convection || !_stage, "for pure diffusion use alternating time steps");
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
        double* tss = elem.time_step_scale();
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

        // apply modal filtering
        if (_use_filter) {
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
              Mat<row_size, Pde::n_update> row_r = Row_rw<Pde::n_update, row_size>::read_row(time_rate[0], ind);
              row_r = filter*row_r;
              Row_rw<Pde::n_update, row_size>::write_row(row_r, time_rate[0], ind, 0.);
            }
          }
        }

        // write update to interior
        for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
          double* to_update = state + ((_compute_residual ? Pde::ref_start : Pde::curr_start) + i_var)*n_qpoint;
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            double det = 1;
            if constexpr (element_t::is_deformed) det = elem_det[i_qpoint];
            to_update[i_qpoint] += _update*time_rate[i_var][i_qpoint]*tss[i_qpoint]/d_pos/det;;
          }
        }
        // *now* we can extrapolate state to faces
        write_face(state, elem.faces);
      }
    }
  };

  /*! \brief Computes the shared numerical flux at the element interfaces.
   * \details Requires that the state has been written to the face storage
   * and replaces the state of both faces with the computed flux.
   */
  template <int n_dim, int row_size>
  class Neighbor : public Kernel<Face_connection<element_t>&>
  {
    using Pde = Pde_templ<n_dim>;
    const Pde eq;
    static constexpr int n_fqpoint = math::pow(row_size, n_dim - 1);

    public:
    template <typename... pde_args>
    Neighbor(pde_args... args) : eq{args...} {}

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

  //! compute the difference between the numerical (average) viscous flux and
  //! the viscous flux on each face in preparation for reconciling the different face fluxes
  template <int n_dim, int row_size>
  class Neighbor_reconcile : public Kernel<Face_connection<element_t>&>
  {
    using Pde = Pde_templ<n_dim>;
    static constexpr int n_fqpoint = math::pow(row_size, n_dim - 1);

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

  //! compute the maximum stable time step
  template <int n_dim, int row_size>
  class Max_dt : public Kernel<element_t&, double>
  {
    using Pde = Pde_templ<n_dim>;
    const Pde eq;
    double max_cfl_c;
    double max_cfl_d;
    Mat<row_size> nodes;
    bool _is_local;

    public:
    template <typename... pde_args>
    Max_dt(const Basis& basis, bool is_local, bool use_filter, double safety_conv, double safety_diff, pde_args... args) :
      eq{args...},
      max_cfl_c{basis.max_cfl()*safety_conv},
      max_cfl_d{-2/basis.min_eig_diffusion()*safety_diff},
      _is_local{is_local}
    {
      for (int i_node = 0; i_node < row_size; ++i_node) nodes(i_node) = basis.node(i_node);
    }

    virtual double operator()(Sequence<element_t&>& elements)
    {
      constexpr int n_qpoint = math::pow(row_size, n_dim);
      // compute the maximum stable time step for all elements and take the minimum
      double dt = std::numeric_limits<double>::max();
      #pragma omp parallel for reduction(min:dt)
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        element_t& elem {elements[i_elem]};
        double* state = elem.stage(0);
        double* art_visc = elem.art_visc_coef();
        double* tss = elem.time_step_scale();
        Mat<math::pow(2, n_dim)> vertex_spacing;
        for (unsigned i_vert = 0; i_vert < vertex_spacing.size(); ++i_vert) {
          vertex_spacing(i_vert) = elem.vertex_time_step_scale(i_vert);
        }
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
        {
          // fetch state variables
          Mat<Pde::n_var> qpoint_state;
          for (int i_var = 0; i_var < Pde::n_var; ++i_var) {
            qpoint_state(i_var) = state[i_var*n_qpoint + i_qpoint];
          }
          // compute speeds
          // get mesh spacing
          Mat<n_dim> coords;
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            coords(i_dim) = nodes((i_qpoint/math::pow(row_size, n_dim - 1 - i_dim))%row_size);
          }
          // compute time step
          double spacing = math::interp(vertex_spacing, coords);
          double scale = 0;
          if constexpr (Pde::has_convection) scale += eq.char_speed(qpoint_state)/max_cfl_c/spacing;
          if constexpr (Pde::is_viscous) scale += eq.diffusivity(qpoint_state, art_visc[i_qpoint])/max_cfl_d/spacing/spacing;
          if (_is_local) tss[i_qpoint] = 1./scale;
          else {
            tss[i_qpoint] = 1.;
            dt = std::min(dt, 1./scale);
          }
        }
      }
      return _is_local ? 1. : dt;
    }
  };
};

}
#endif
