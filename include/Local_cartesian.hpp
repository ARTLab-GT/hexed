#ifndef HEXED_LOCAL_CARTESIAN_HPP_
#define HEXED_LOCAL_CARTESIAN_HPP_

#include "Vector_view.hpp"
#include "Element.hpp"
#include "Basis.hpp"
#include "kernel_factory.hpp"
#include "math.hpp"
#include "Derivative.hpp"
#include "Write_face.hpp"
#include "thermo.hpp"

namespace hexed
{

class Nd_index
{
  int i_outer = 0;
  int i_inner = 0;
  int i_fq = 0;
  public:
  const int n_dim;
  const int row_size;
  const int i_dim;
  const int n_qpoint;
  const int n_fqpoint;
  const int stride;
  constexpr Nd_index(int nd, int rs, int id)
  : n_dim{nd}, row_size{rs}, i_dim{id},
    n_qpoint{custom_math::pow(row_size, n_dim)},
    n_fqpoint{n_qpoint/row_size},
    stride{custom_math::pow(row_size, n_dim - 1 - i_dim)}
  {}
  constexpr void operator++()
  {
    ++i_fq;
    ++i_inner;
    if (i_inner == stride) {
      i_inner = 0;
      ++i_outer;
    }
  }
  constexpr operator bool() const {return i_fq < n_fqpoint;}
  constexpr int i_face_qpoint() const {return i_fq;}
  constexpr int i_qpoint(int i_node) const {return i_outer*stride*row_size + i_inner + i_node*stride;}
};

template <int rows, int cols = 1>
using Mat = Eigen::Matrix<double, rows, cols>;

template <int rows = Eigen::Dynamic, int cols = Eigen::Dynamic>
using Map = Eigen::Map<Mat<rows, cols>>;

template <int n_var, int row_size>
class Numerics
{
  public:
  typedef Mat<row_size, n_var> Row;
  typedef Mat<2, n_var> Bound;

  Numerics() = delete;

  static Row read_row(const double* data, Nd_index ind)
  {
    Row r;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_row = 0; i_row < row_size; ++i_row) {
        r(i_row, i_var) = data[i_var*ind.n_qpoint + ind.i_qpoint(i_row)];
      }
    }
    return r;
  }

  static void write_row(Row w, double* data, Nd_index ind, double coef)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_row = 0; i_row < row_size; ++i_row) {
        double& d = data[i_var*ind.n_qpoint + ind.i_qpoint(i_row)];
        d = coef*d + w(i_row, i_var);
      }
    }
  }

  static Bound read_bound(const double* face, Nd_index ind)
  {
    Bound b;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        b(is_positive, i_var) = face[((ind.i_dim*2 + is_positive)*n_var + i_var)*ind.n_fqpoint + ind.i_face_qpoint()];
      }
    }
    return b;
  }

  static void write_bound(Bound b, double* face, Nd_index ind)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        face[((ind.i_dim*2 + is_positive)*n_var + i_var)*ind.n_fqpoint + ind.i_face_qpoint()] = b(is_positive, i_var);
      }
    }
  }
};

class Row_ops
{
  public:
  const int n_dim;
  const int row_size;
  const int n_qpoint;
  const int n_fqpoint;
  constexpr int stride(int i_dim) const {return custom_math::pow(row_size, n_dim - 1 - i_dim);}
  constexpr Row_ops(int nd, int rs) :
    n_dim{nd},
    row_size{rs},
    n_qpoint{custom_math::pow(row_size, n_dim)},
    n_fqpoint{n_qpoint/row_size}
  {}

  auto row(double* data, int n_var, int i_fqpoint, int i_dim) const
  {
    Map<> m(data, n_qpoint, n_var);
    int s = stride(i_dim);
    return m(Eigen::seqN(i_fqpoint%s + row_size*s*(i_fqpoint/s), row_size, s), Eigen::all);
  }

  auto boundary(double* face, int n_var, int i_fqpoint, int i_dim, int offset = 0, int n_var_read = 0) const
  {
    Map<> m(face, n_var*n_fqpoint, 2*n_dim);
    Eigen::MatrixXd m1 = m(Eigen::seqN(i_fqpoint + offset*n_fqpoint, (n_var_read ? n_var_read : n_var), n_fqpoint), Eigen::seqN(2*i_dim, 2));
    Eigen::MatrixXd m2 = m1.transpose();
    return m2;
  }
};

#define HEXED_ASSERT_THERM_ADMIS \
  HEXED_ASSERT(state(n_dim) > 0, "nonpositive density"); \
  HEXED_ASSERT(state(n_dim + 1) >= 0, "negative energy"); \
  HEXED_ASSERT(pres >= 0, "negative pressure"); \

template <int n_dim>
class Physics
{
  public:
  Physics() = delete;
  static constexpr int n_var = n_dim + 2;
  static constexpr double heat_rat = 1.4;

  static constexpr double pressure(Mat<n_var> state)
  {
    double mmtm_sq = 0.;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      mmtm_sq += (state(j_dim))*(state(j_dim));
    }
    return (heat_rat - 1.)*((state(n_dim + 1)) - 0.5*mmtm_sq/(state(n_dim)));
  }

  static constexpr Mat<n_var> flux(Mat<n_var> state, Mat<n_dim> normal, int i_dim)
  {
    Mat<n_var> f;
    f(n_dim) = 0.;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      f(n_dim) += state(j_dim)*normal(j_dim);
    }
    double scaled = f(n_dim)/state(n_dim);
    double pres = pressure(state);
    HEXED_ASSERT_THERM_ADMIS
    f(n_var - 1) = (state(n_dim + 1) + pres)*scaled;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      f(j_dim) = state(j_dim)*scaled + pres*normal(j_dim);
    }
    return f;
  }
};

template <typename element_t>
class Transform_data
{
  public:
  double* normal;
  double* jac_det;
  Transform_data(element_t&, double*);
};

template<>
Transform_data<Element>::Transform_data(Element& elem, double* dummy) :
  normal{dummy},
  jac_det{dummy}
{}

template<>
Transform_data<Deformed_element>::Transform_data(Deformed_element& elem, double* dummy) :
  normal{elem.reference_level_normals()},
  jac_det{elem.jacobian_determinant()}
{}

template <typename element_t>
class Spatial
{
  public:
  Spatial() = delete;

  template <int n_dim, int row_size>
  class Local : public Kernel<element_t&>
  {
    using Phys = Physics<n_dim>;
    using Num = Numerics<Phys::n_var, row_size>;
    static constexpr Row_ops ro{n_dim, row_size};
    Derivative<row_size> derivative;
    Write_face<n_dim, row_size> write_face;
    double update;
    double curr;
    double ref;
    const double heat_rat;
    static constexpr int n_qpoint = custom_math::pow(row_size, n_dim);
    static constexpr int n_face_dof = Phys::n_var*n_qpoint/row_size;
    double dummy [n_dim][n_dim][n_qpoint];

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
    {
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            dummy[i_dim][j_dim][i_qpoint] = (i_dim == j_dim);
          }
        }
      }
    }

    virtual void operator()(Sequence<element_t&>& elements)
    {
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        double* active_state = elem.stage(0);
        double* reference_state = active_state + Phys::n_var*n_qpoint;
        double* face = elem.face();
        double* tss = elem.time_step_scale();
        double d_pos = elem.nominal_size();
        Transform_data dat(elem, dummy[0][0]);
        double time_rate [Phys::n_var][n_qpoint] {};

        // compute update
        for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
          for (int i_fqpoint = 0; i_fqpoint < ro.n_fqpoint; ++i_fqpoint) {
            Mat<row_size, Phys::n_var> row_r = ro.row(active_state, Phys::n_var, i_fqpoint, i_dim);
            Mat<row_size, n_dim      > row_n = ro.row(dat.normal + i_dim*n_dim*n_qpoint, n_dim, i_fqpoint, i_dim);
            Mat<row_size, Phys::n_var> flux;
            for (int i_row = 0; i_row < row_size; ++i_row) {
              flux(i_row, Eigen::all) = Phys::flux(row_r(i_row, Eigen::all), row_n(i_row, Eigen::all), i_dim);
            }
            Mat<2, Phys::n_var> face_r = ro.boundary(face, Phys::n_var, i_fqpoint, i_dim);
            ro.row(time_rate[0], Phys::n_var, i_fqpoint, i_dim) -= derivative(flux, face_r);
          }
        }

        // write the updated solution
        for (int i_var = 0; i_var < Phys::n_var; ++i_var) {
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
            const int i_dof = i_var*n_qpoint + i_qpoint;
            active_state[i_dof] = update*time_rate[i_var][i_qpoint]/d_pos*tss[i_qpoint]/dat.jac_det[i_qpoint]
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
