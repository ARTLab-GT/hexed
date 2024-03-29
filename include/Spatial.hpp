#ifndef HEXED_SPATIAL_HPP_
#define HEXED_SPATIAL_HPP_

#include "Sequence.hpp"
#include "Kernel_element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Row_rw.hpp"
#include "Face_permutation.hpp"
#include "Refined_face.hpp"

namespace hexed
{

//! class to contain all the kernels in a scope
//! parameterized by the type of element (deformed/cartesian) and the PDE
template <template<int, int> typename Pde_templ, bool is_deformed>
class Spatial
{
  public:
  //! this class has no non-static data members, so there is no reason to construct it
  Spatial() = delete;
  Spatial(const Spatial&) = delete;
  Spatial(Spatial&&) = delete;

  //! extrapolates the values in the interior of an element to the faces
  template <int n_dim, int row_size>
  class Write_face : public Kernel<Kernel_element&>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    Pde _eq;
    const Eigen::Matrix<double, 2, row_size> boundary;

    public:
    template <typename... pde_args>
    Write_face(const Basis& basis, pde_args... args) : _eq(args...), boundary{basis.boundary()} {}

    //! apply to a single element
    void operator()(const double* read, std::array<double*, 6> faces)
    {
      constexpr int n_qpoint = math::pow(row_size, n_dim);
      double extrap [Pde::n_extrap][n_qpoint];
      // fetch the extrapolation variables and store them in `time_rate`
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        Mat<Pde::n_extrap> grad_vars = _eq.fetch_extrap(n_qpoint, read + i_qpoint);
        for (int i_var = 0; i_var < Pde::n_extrap; ++i_var) extrap[i_var][i_qpoint] = grad_vars(i_var);
      }
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
          auto row_r = Row_rw<Pde::n_extrap, row_size>::read_row(extrap[0], ind);
          Mat<2, Pde::n_extrap> bound = boundary*row_r;
          Row_rw<Pde::n_extrap, row_size>::write_bound(bound, faces, ind);
        }
      }
    }

    //! apply to a sequence of elements
    virtual void operator()(Sequence<Kernel_element&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
        Kernel_element& elem {elements[i_elem]};
        double* read = elem.state();
        std::array<double*, 6> faces {};
        for (int i_face = 0; i_face < 2*n_dim; ++i_face) faces[i_face] = elem.face(i_face, false);
        operator()(read, faces);
      }
    }
  };

  template<int n_dim, int row_size>
  class Face_permutation : public Face_permutation_dynamic
  {
    static constexpr int n_qpoint = math::pow(row_size, n_dim - 1);
    static constexpr int n_var = Pde_templ<n_dim, row_size>::n_extrap;
    Connection_direction dir;
    double* tgt;

    /*
     * For 3d connections, the dimensions corresponding to major/minor indices may not match.
     * To solve this problem, this function transposes the major/minor axes if necessary.
     */
    virtual void transpose()
    {
      if (dir.transpose()) {
        if constexpr (n_dim == 3) {
          for (int i_var = 0; i_var < n_var; ++i_var) {
            Eigen::Map<Eigen::Matrix<double, row_size, row_size>> rows {tgt + i_var*n_qpoint};
            rows.transposeInPlace();
          }
        }
      }
    }

    /* For 2d and 3d connections, the ordering of quadrature points in the plane of the
     * normal vectors might be opposite.
     * To solve this problem, this function reverses the order in `i_dim[0]` if necessary.
     */
    virtual void flip()
    {
      if (dir.flip_tangential()) {
        if constexpr (n_dim == 3) {
          // figure out if the dimension we want to reverse is the major or minor axis of face 1
          bool colwise = (dir.i_dim[0] > 3 - dir.i_dim[0] - dir.i_dim[1]) != dir.transpose();
          for (int i_var = 0; i_var < n_var; ++i_var) {
            Eigen::Map<Eigen::Matrix<double, row_size, row_size>> rows {tgt + i_var*n_qpoint};
            if (colwise) rows.colwise().reverseInPlace();
            else         rows.rowwise().reverseInPlace();
          }
        } else if constexpr (n_dim == 2) {
          // for 2d, there is no major/minor axis distinction
          Eigen::Map<Eigen::Matrix<double, row_size, n_var*n_qpoint/row_size>> rows {tgt};
          rows.colwise().reverseInPlace();
        }
      }
    }

    public:
    // `direction` defines the manner in which the faces are connected
    // `target` points to the data for face 1
    // note that this might be the actual `Element::face()` data, or some temporary storage allocated by the caller,
    // but it should correspond to face 1
    Face_permutation(Connection_direction direction, double* target)
    : dir{direction}, tgt{target} {}
    // reordering can consist of order reversal and/or transpose operations
    virtual void match_faces() {transpose(); flip();}
    virtual void restore()     {flip(); transpose();}
    typedef std::unique_ptr<Face_permutation_dynamic> ptr_t;
  };

  /*
   * See `Refined_face.hpp`.
   * Transforms data from coarse face to fine faces by polynomial interpolation.
   * This interpolation is exact (and therefore conservative).
   */
  template <int n_dim, int row_size>
  class Prolong_refined : public Kernel<Refined_face&>
  {
    static constexpr int n_var = Pde_templ<n_dim, row_size>::n_extrap;
    const Eigen::Matrix<double, row_size, row_size> prolong_mat [2];
    bool scl;
    bool off;

    public:
    Prolong_refined(const Basis& basis, bool scale = false, bool offset = false) :
      prolong_mat{basis.prolong(0), basis.prolong(1)},
      scl{scale},
      off{offset}
    {}

    virtual void operator()(Sequence<Refined_face&>& ref_faces)
    {
      constexpr int n_face = math::pow(2, n_dim - 1);
      constexpr int nfq = math::pow(row_size, n_dim - 1);

      #pragma omp parallel for
      for (int i_ref_face = 0; i_ref_face < ref_faces.size(); ++i_ref_face)
      {
        auto& ref_face {ref_faces[i_ref_face]};
        double* coarse {ref_face.coarse + off*(n_dim + 2)*nfq};
        const auto str = ref_face.stretch;
        // update number of faces to reflect any face stretching
        int nf = n_face;
        for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim) nf /= 1 + str[i_dim];
        for (int i_face = 0; i_face < nf; ++i_face)
        {
          double* fine {ref_face.fine[i_face] + off*(n_dim + 2)*nfq};
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            double* var_face {fine + i_var*nfq};
            // initialize fine face to be equal to coarse face
            for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
              var_face[i_qpoint] = coarse[i_var*nfq + i_qpoint];
            }
            // interpolate one dimension at a time, in-place
            for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim)
            {
              if (str[i_dim]) for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) var_face[i_qpoint] *= 1 + scl;
              else
              {
                const int pow {n_dim - 2 - i_dim};
                const int face_stride {str[n_dim - 2] ? 1 : math::pow(2, pow)};
                const int qpoint_stride {math::pow(row_size, pow)};
                const int i_half {(i_face/face_stride)%2}; // is this face covering the upper or lower half of the coarse face with respect to the current dimension?
                for (int i_outer = 0; i_outer < nfq/(row_size*qpoint_stride); ++i_outer)
                {
                  for (int i_inner = 0; i_inner < qpoint_stride; ++i_inner)
                  {
                    Eigen::Matrix<double, row_size, 1> row;
                    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                      row(i_qpoint) = var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner];
                    }
                    row = prolong_mat[i_half]*row;
                    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                      var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner] = row(i_qpoint);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  };

  /*
   * See `Refined_face.hpp`.
   * Transforms data from fine faces to coarse face by orthogonal projection.
   * This projection is conservative (in the sense that the integral of the face data does not change)
   * but of course not always exact.
   */
  template <int n_dim, int row_size>
  class Restrict_refined : public Kernel<Refined_face&>
  {
    static constexpr int n_var = Pde_templ<n_dim, row_size>::n_extrap;
    const Eigen::Matrix<double, row_size, row_size> restrict_mat [2];
    bool scl;
    bool off;

    public:
    Restrict_refined(const Basis& basis, bool scale = true, bool offset = false) :
      restrict_mat{basis.restrict(0), basis.restrict(1)},
      scl{scale},
      off{offset}
    {}

    virtual void operator()(Sequence<Refined_face&>& ref_faces)
    {
      constexpr int n_face = math::pow(2, n_dim - 1);
      constexpr int nfq = math::pow(row_size, n_dim - 1);

      #pragma omp parallel for
      for (int i_ref_face = 0; i_ref_face < ref_faces.size(); ++i_ref_face)
      {
        auto& ref_face {ref_faces[i_ref_face]};
        double* coarse {ref_face.coarse + off*(n_dim + 2)*nfq};
        for (int i_dof = 0; i_dof < n_var*nfq; ++i_dof) coarse[i_dof] = 0.;
        auto str = ref_face.stretch;
        // update number of faces to reflect any face stretching
        int nf = n_face;
        for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim) nf /= 1 + str[i_dim];
        for (int i_face = 0; i_face < nf; ++i_face)
        {
          double* fine {ref_face.fine[i_face] + off*(n_dim + 2)*nfq};
          for (int i_var = 0; i_var < n_var; ++i_var)
          {
            double* var_face {fine + i_var*nfq};
            for (int i_dim = 0; i_dim < n_dim - 1; ++i_dim)
            {
              if (str[i_dim])
              {
                for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) var_face[i_qpoint] /= 1 + scl;
              }
              else
              {
                const int pow {n_dim - 2 - i_dim};
                const int face_stride {str[n_dim - 2] ? 1 : math::pow(2, pow)};
                const int qpoint_stride {math::pow(row_size, pow)};
                const int i_half {(i_face/face_stride)%2}; // is this face covering the upper or lower half of the coarse face with respect to the current dimension?
                for (int i_outer = 0; i_outer < nfq/(row_size*qpoint_stride); ++i_outer)
                {
                  for (int i_inner = 0; i_inner < qpoint_stride; ++i_inner)
                  {
                    Eigen::Matrix<double, row_size, 1> row;
                    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                      row(i_qpoint) = var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner];
                    }
                    row = restrict_mat[i_half]*row;
                    for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
                      var_face[(i_outer*row_size + i_qpoint)*qpoint_stride + i_inner] = row(i_qpoint);
                    }
                  }
                }
              }
            }
            for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
              coarse[i_var*nfq + i_qpoint] += var_face[i_qpoint];
            }
          }
        }
      }
    }
  };

  /*! \brief Performs the update to the element state after the shared numerical flux has been computed.
   * \details In other words, performs the part of the algorithm that only uses local data
   * (or in other words does not directly depend on neighboring elements).
   * Note that face data is updated to reflect the updated interior state.
   */
  template <int n_dim, int row_size>
  class Local : public Kernel<Kernel_element&>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    const Pde _eq;
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
      _eq(args...),
      derivative{basis},
      boundary{basis.boundary()},
      write_face(basis, args...),
      filter{basis.filter()},
      _update{stage ? dt*basis.step_ratio() : dt},
      _stage{stage},
      _compute_residual{compute_residual},
      _use_filter{use_filter}
    {
      HEXED_ASSERT(!(Pde::has_diffusion & _stage), "two-stage stabilization is not applicable to diffusion equations");
      HEXED_ASSERT(!(_stage && _compute_residual), "residual calculation is a single-stage operation");
    }

    virtual void operator()(Sequence<Kernel_element&>& elements)
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
        double* state = elem.state();
        std::array<double*, 6> faces;
        for (int i_face = 0; i_face < 2*n_dim; ++i_face) faces[i_face] = elem.face(i_face, false);
        std::array<double*, 6> visc_faces;
        for (int i_face = 0; i_face < 2*n_dim; ++i_face) visc_faces[i_face] = elem.face(i_face, true);
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [2][std::max(Pde::n_update, Pde::n_extrap - Pde::n_extrap/2)][n_qpoint] {}; // first part contains convective time derivative, second part diffusive
        // only need the next 2 for deformed elements
        double* nrml = nullptr; // reference level normals
        double* elem_det = nullptr; // jacobian determinant
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        std::array<double*, 6> face_nrml;
        double visc_storage [n_dim][Pde::n_extrap][n_qpoint] {}; // for viscous PDEs, will be used to store gradients and later viscous flux
        #pragma GCC diagnostic pop
        if constexpr (is_deformed) {
          elem_det = elem.jacobian_determinant();
          nrml = elem.reference_level_normals();
          for (int i_face = 0; i_face < 2*n_dim; ++i_face) {
            face_nrml[i_face] = elem.kernel_face_normal(i_face);
            if (!face_nrml[i_face]) face_nrml[i_face] = cartesian_normal[i_face/2][0];
          }
        }

        // compute gradient (times jacobian determinant, cause that's easier)
        if constexpr (Pde::has_diffusion)
        {
          static_assert(Pde::n_extrap >= Pde::n_update);
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            Mat<Pde::n_extrap> grad_vars = _eq.fetch_extrap(n_qpoint, state + i_qpoint);
            for (int i_var = 0; i_var < Pde::n_extrap; ++i_var) (&time_rate[0][0][i_qpoint])[i_var*n_qpoint] = grad_vars(i_var);
          }
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
              // fetch state data
              auto row_r = Row_rw<Pde::n_extrap, row_size>::read_row(time_rate[0][0], ind);
              auto face_state = Row_rw<Pde::n_extrap, row_size>::read_bound(visc_faces, ind);
              if constexpr (is_deformed) {
                // fetch normal data
                auto row_n = Row_rw<n_dim, row_size>::read_row(nrml + i_dim*n_dim*n_qpoint, ind);
                auto face_n = Row_rw<n_dim, row_size>::read_bound(face_nrml, ind);
                // differentiate and write to temporary storage
                for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
                  for (int i_var = 0; i_var < Pde::n_extrap; ++i_var) {
                    Mat<row_size, 1> row = row_n(Eigen::all, j_dim).cwiseProduct(row_r(Eigen::all, i_var));
                    Mat<2, 1> bound = face_n(Eigen::all, j_dim).cwiseProduct(face_state(Eigen::all, i_var));
                    Row_rw<1, row_size>::write_row(derivative(row, bound)/d_pos, visc_storage[j_dim][i_var], ind, 1.);
                  }
                }
              } else {
                // differentiate and write to temporary storage
                Row_rw<Pde::n_extrap, row_size>::write_row(derivative(row_r, face_state)/d_pos, visc_storage[i_dim][0], ind, 0);
              }
            }
          }
          for (int i = 0; i < 2*Pde::n_update*n_qpoint; ++i) (&time_rate[0][0][0])[i] = 0;
        }

        // compute flux
        double flux [n_dim][Pde::n_update][n_qpoint];
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          typename Pde::Computation<n_dim> comp(_eq);
          comp.fetch_state(n_qpoint, state + i_qpoint);
          if constexpr (is_deformed) {
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int j_dim = 0; j_dim < n_dim; ++j_dim) comp.normal(j_dim, i_dim) = nrml[(i_dim*n_dim + j_dim)*n_qpoint + i_qpoint];
            }
          }
          if constexpr (Pde::has_convection) {
            comp.compute_flux_conv();
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int i_var = 0; i_var < Pde::n_update; ++i_var) flux[i_dim][i_var][i_qpoint] = comp.flux_conv(i_var, i_dim);
            }
          }
          if constexpr (Pde::has_diffusion) {
            // fetch gradient
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int i_var = 0; i_var < Pde::n_extrap; ++i_var) {
                comp.gradient(i_var, i_dim) = visc_storage[i_dim][i_var][i_qpoint];
              }
            }
            if constexpr (is_deformed) comp.gradient /= elem_det[i_qpoint]; // divide by the determinant to get the actual gradient (see above)
            // compute flux and write to temporary storage
            comp.compute_flux_diff();
            for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
              for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
                visc_storage[i_dim][i_var][i_qpoint] = comp.flux_diff(i_var, i_dim);
              }
            }
          }
          if constexpr (Pde::has_source) if (!_stage) {
            comp.compute_source();
            double mult = d_pos;
            if constexpr (is_deformed) mult *= elem_det[i_qpoint];
            for (int i_var = 0; i_var < Pde::n_update; ++i_var) time_rate[1][i_var][i_qpoint] = mult*comp.source(i_var);
          }
        }

        // compute residual
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            Mat<row_size, Pde::n_update> row_f;
            Mat<2, Pde::n_update> face_f;
            if constexpr (Pde::has_convection) {
              // fetch row data
              row_f = Row_rw<Pde::n_update, row_size>::read_row(flux[i_dim][0], ind);
              // fetch face data
              face_f = Row_rw<Pde::n_update, row_size>::read_bound(faces, ind);
              // differentiate and write to temporary storage
              Row_rw<Pde::n_update, row_size>::write_row(-derivative(row_f, face_f), time_rate[0][0], ind, 1.);
            }
            // compute viscous update
            if constexpr (Pde::has_diffusion) {
              row_f = Row_rw<Pde::n_update, row_size>::read_row(visc_storage[i_dim][0], ind);
              face_f = boundary*row_f;
              // write viscous row_f to faces to enable calculation of the numerical row_f
              for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
                for (int is_positive : {0, 1}) {
                  visc_faces[ind.i_dim*2 + is_positive][i_var*ind.n_fqpoint + ind.i_face_qpoint()] = face_f(is_positive, i_var);
                }
              }
              Row_rw<Pde::n_update, row_size>::write_row(-derivative(row_f), time_rate[1][0], ind, 1.);
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
        double* ref_state = elem.residual_cache();
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          Mat<Pde::n_update> update;
          update.setZero();
          double mult = _update*tss[i_qpoint]/d_pos;
          if constexpr (is_deformed) mult /= elem_det[i_qpoint];
          for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
            double u = time_rate[0][i_var][i_qpoint];
            if (_stage) {
              u -= ref_state[i_var*n_qpoint + i_qpoint];
            } else {
              if constexpr (Pde::has_convection) ref_state[i_var*n_qpoint + i_qpoint] = u;
              if constexpr (Pde::has_diffusion || Pde::has_source) u += time_rate[1][i_var][i_qpoint];
            }
            u *= mult;
            if (_compute_residual) ref_state[i_var*n_qpoint + i_qpoint] = u;
            else update(i_var) = u;
          }
          _eq.write_update(update, n_qpoint, state + i_qpoint, !Pde::has_diffusion && !_stage);
        }

        // write updated state to face storage.
        // For viscous, don't bother since we still have to add the numerical flux term
        if constexpr (!Pde::has_diffusion) write_face(state, faces);
      }
    }
  };

  //! account for the difference between the real and numerical viscous fluxes to enforce conservation
  template <int n_dim, int row_size>
  class Reconcile_ldg_flux : public Kernel<Kernel_element&>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    const Pde _eq;
    static constexpr int n_qpoint = math::pow(row_size, n_dim);
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    Mat<row_size, row_size> filter;
    double _update;
    int _stage;
    bool _compute_residual;
    bool _use_filter;

    public:
    template <typename... pde_args>
    Reconcile_ldg_flux(const Basis& basis, double dt, int which_stage, bool compute_residual, bool use_filter, pde_args... args) :
      _eq(args...),
      derivative{basis},
      write_face(basis, args...),
      filter{basis.filter()},
      _update{dt},
      _stage{which_stage},
      _compute_residual{compute_residual},
      _use_filter{use_filter}
    {
      HEXED_ASSERT(!(Pde::has_diffusion & _stage), "two-stage stabilization is not applicable to diffusion equations");
      HEXED_ASSERT(Pde::has_convection || !_stage, "for pure diffusion use alternating time steps");
    }

    virtual void operator()(Sequence<Kernel_element&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* state = elem.state();
        std::array<double*, 6> visc_faces;
        for (int i_face = 0; i_face < 2*n_dim; ++i_face) visc_faces[i_face] = elem.face(i_face, true);
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        double time_rate [Pde::n_update][n_qpoint] {};
        double* elem_det = nullptr;
        if constexpr (is_deformed) {
          elem_det = elem.jacobian_determinant();
        }
        // fetch the flux difference and lift it to the interior temporary storage
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (Row_index ind(n_dim, row_size, i_dim); ind; ++ind) {
            auto bound_f = Row_rw<Pde::n_update, row_size>::read_bound(visc_faces, ind);
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
        double* to_update = _compute_residual ? elem.residual_cache() : state;
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          Mat<Pde::n_update> update;
          double mult = _update*tss[i_qpoint]/d_pos;
          if constexpr (is_deformed) mult /= elem_det[i_qpoint];
          for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
            update(i_var) = time_rate[i_var][i_qpoint]*mult;
          }
          _eq.write_update(update, n_qpoint, to_update + i_qpoint, true);
        }
        // *now* we can extrapolate state to faces
        std::array<double*, 6> faces;
        for (int i_face = 0; i_face < 2*n_dim; ++i_face) faces[i_face] = elem.face(i_face, false);
        write_face(state, faces);
      }
    }
  };

  /*! \brief Computes the shared numerical flux at the element interfaces.
   * \details Requires that the state has been written to the face storage
   * and replaces the state of both faces with the computed flux.
   */
  template <int n_dim, int row_size>
  class Neighbor : public Kernel<Kernel_connection&>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    const Pde _eq;
    static constexpr int n_fqpoint = math::pow(row_size, n_dim - 1);
    const int _stage;

    public:
    template <typename... pde_args>
    Neighbor(int i_stage, pde_args... args) : _eq(args...), _stage{i_stage} {}

    virtual void operator()(Sequence<Kernel_connection&>& connections)
    {
      #pragma omp parallel for
      for (int i_con = 0; i_con < connections.size(); ++i_con)
      {
        auto& con = connections[i_con];
        auto dir = con.get_direction();
        double face [2 + 2*Pde::has_diffusion][Pde::n_extrap*n_fqpoint] {}; // copying face data to temporary stack storage improves efficiency
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wunused-but-set-variable"
        double face_nrml [n_dim*n_fqpoint]; // only set for deformed
        #pragma GCC diagnostic pop
        int sign [2] {1, 1}; // records whether the normal vector on each side needs to be flipped to obey sign convention
        // fetch face data
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state(i_side, false);
          for (int i_dof = 0; i_dof < Pde::n_extrap*n_fqpoint; ++i_dof) {
            face[i_side][i_dof] = f[i_dof];
          }
        }
        Face_permutation<n_dim, row_size> perm(dir, face[1]); // only used for deformed
        if constexpr (is_deformed) {
          perm.match_faces(); // if order of quadrature points on both faces does not match, reorder face 1 to match face 0
          for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
          double* n = con.normal();
          for (int i_dof = 0; i_dof < n_dim*n_fqpoint; ++i_dof) {
            face_nrml[i_dof] = n[i_dof];
          }
        }
        // compute flux
        for (int i_qpoint = 0; i_qpoint < n_fqpoint; ++i_qpoint)
        {
          // compute average face state for LDG scheme
          if constexpr (Pde::has_diffusion) {
            for (int i_var = 0; i_var < Pde::n_extrap; ++i_var) {
              double avg = .5*(face[0][i_var*n_fqpoint + i_qpoint] + face[1][i_var*n_fqpoint + i_qpoint]);
              for (int i_side = 0; i_side < 2; ++i_side) {
                face[2 + i_side][i_var*n_fqpoint + i_qpoint] = avg;
              }
            }
          }
          if constexpr (Pde::has_convection) {
            // fetch data
            typename Pde::Computation<1> comp [2] {_eq, _eq};
            if constexpr (is_deformed) {
              for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
                comp[0].normal(i_dim) = sign[0]*face_nrml[i_dim*n_fqpoint + i_qpoint];
              }
            } else comp[0].normal.setUnit(dir.i_dim[0]);
            comp[1].normal = comp[0].normal;
            for (int i_side = 0; i_side < 2; ++i_side) {
              comp[i_side].fetch_extrap_state(n_fqpoint, face[i_side] + i_qpoint);
            }
            // compute flux
            for (int i_side = 0; i_side < 2; ++i_side) {
              comp[i_side].compute_flux_conv();
              comp[i_side].compute_char_speed();
            }
            Mat<Pde::n_update> flux = .5*(
              comp[0].flux_conv + comp[1].flux_conv
              + std::max(comp[0].char_speed, comp[1].char_speed)*comp[0].normal.norm()*(comp[0].update_state - comp[1].update_state)
            );
            for (int i_side = 0; i_side < 2; ++i_side) {
              for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
                face[i_side][i_var*n_fqpoint + i_qpoint] = sign[i_side]*flux(i_var);
              }
            }
          }
        }
        if constexpr (is_deformed) {
          perm.restore(); // restore data of face 1 to original order
          if constexpr (Pde::has_diffusion) Face_permutation<n_dim, row_size>(dir, face[3]).restore();
        }
        // write data to actual face storage on heap
        for (int i_side = 0; i_side < 2; ++i_side) {
          // write flux
          if constexpr (Pde::has_convection) {
            double* f = con.state(i_side, false);
            for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
              f[i_dof] = face[i_side][i_dof];
            }
          }
          // write average face state
          if constexpr (Pde::has_diffusion) {
            double* f = con.state(i_side, true);
            for (int i_dof = 0; i_dof < Pde::n_extrap*n_fqpoint; ++i_dof) {
              f[i_dof] = face[2 + i_side][i_dof];
            }
          }
        }
      }
    }
  };

  //! compute the difference between the numerical (average) viscous flux and
  //! the viscous flux on each face in preparation for reconciling the different face fluxes
  template <int n_dim, int row_size>
  class Neighbor_reconcile : public Kernel<Kernel_connection&>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    static constexpr int n_fqpoint = math::pow(row_size, n_dim - 1);

    public:
    virtual void operator()(Sequence<Kernel_connection&>& connections)
    {
      #pragma omp parallel for
      for (int i_con = 0; i_con < connections.size(); ++i_con)
      {
        auto& con = connections[i_con];
        auto dir = con.get_direction();
        double face [2][(n_dim + 2)*n_fqpoint]; // copying face data to temporary stack storage improves efficiency
        int sign [2] {1, 1}; // records whether the normal vector on each side needs to be flipped to obey sign convention
        // fetch face data
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state(i_side, true);
          for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
            face[i_side][i_dof] = f[i_dof];
          }
        }
        Face_permutation<n_dim, row_size> perm(dir, face[1]); // only used for deformed
        if constexpr (is_deformed) {
          perm.match_faces(); // if order of quadrature points on both faces does not match, reorder face 1 to match face 0
          for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
        }
        // compute flux
        for (int i_qpoint = 0; i_qpoint < n_fqpoint; ++i_qpoint) {
          for (int i_var = 0; i_var < Pde::n_update; ++i_var) {
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
        if constexpr (is_deformed) perm.restore(); // restore data of face 1 to original order
        // write data to actual face storage on heap
        for (int i_side = 0; i_side < 2; ++i_side) {
          double* f = con.state(i_side, true);
          for (int i_dof = 0; i_dof < Pde::n_update*n_fqpoint; ++i_dof) {
            f[i_dof] = face[i_side][i_dof];
          }
        }
      }
    }
  };

  //! compute the maximum stable time step
  template <int n_dim, int row_size>
  class Max_dt : public Kernel<Kernel_element&, double>
  {
    using Pde = Pde_templ<n_dim, row_size>;
    const Pde _eq;
    double max_cfl_c;
    double max_cfl_d;
    Mat<row_size> nodes;
    bool _is_local;

    public:
    template <typename... pde_args>
    Max_dt(const Basis& basis, bool is_local, bool use_filter, double safety_conv, double safety_diff, pde_args... args) :
      _eq(args...),
      max_cfl_c{basis.max_cfl()*safety_conv},
      max_cfl_d{-2/basis.min_eig_diffusion()*safety_diff},
      _is_local{is_local}
    {
      for (int i_node = 0; i_node < row_size; ++i_node) nodes(i_node) = basis.node(i_node);
    }

    virtual double operator()(Sequence<Kernel_element&>& elements)
    {
      constexpr int n_qpoint = math::pow(row_size, n_dim);
      // compute the maximum stable time step for all elements and take the minimum
      double dt = std::numeric_limits<double>::max();
      #pragma omp parallel for reduction(min:dt)
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        Kernel_element& elem {elements[i_elem]};
        double* state = elem.state();
        double* tss = elem.time_step_scale();
        Mat<math::pow(2, n_dim)> vertex_spacing;
        for (unsigned i_vert = 0; i_vert < vertex_spacing.size(); ++i_vert) {
          vertex_spacing(i_vert) = elem.vertex_time_step_scale(i_vert);
        }
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
        {
          // get mesh spacing
          Mat<n_dim> coords;
          for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
            coords(i_dim) = nodes((i_qpoint/math::pow(row_size, n_dim - 1 - i_dim))%row_size);
          }
          double spacing = math::interp(vertex_spacing, coords);
          // fetch state
          typename Pde::Computation<n_dim> comp(_eq);
          comp.fetch_state(n_qpoint, state + i_qpoint);
          // compute time step
          double scale = 0;
          if constexpr (Pde::has_convection) {
            comp.compute_char_speed();
            scale += comp.char_speed/max_cfl_c/spacing;
          }
          if constexpr (Pde::has_diffusion) {
            comp.compute_diffusivity();
            scale += comp.diffusivity/max_cfl_d/spacing/spacing;
          }
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
