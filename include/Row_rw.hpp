#ifndef HEXED_ROW_RW_HPP_
#define HEXED_ROW_RW_HPP_

#include "Row_index.hpp"

namespace hexed
{

/*
 * "row read/write".
 * Can read/write to/from a row of quadrature points in an element
 * and/or the associated face quadrature points.
 * This class has only static members, so it is basically just a templated namespace.
 * Note that the dimensionality of the element data is determined
 * by the `Row_index` argument to the member functions, not inferred from `n_var`
 * (which is arbitrary regardless of dimensionality).
 */
template <int n_var, int row_size>
class Row_rw
{
  public:
  typedef Mat<row_size, n_var> Row;
  typedef Mat<2, n_var> Bound;

  // since this class is used as a templated namespace, its constructors are deleted
  Row_rw() = delete;
  Row_rw(const Row_rw&) = delete;
  Row_rw(Row_rw&&) = delete;

  // reads from a row of quadrature points
  static Row read_row(const double* data, Row_index ind)
  {
    Row r;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_row = 0; i_row < row_size; ++i_row) {
        r(i_row, i_var) = data[i_var*ind.n_qpoint + ind.i_qpoint(i_row)];
      }
    }
    return r;
  }

  // multiples a row of values by `coef` and then adds the values in `w` to it
  // (basically a row-wise axpy)
  static void write_row(Row w, double* data, Row_index ind, double coef)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_row = 0; i_row < row_size; ++i_row) {
        double& d = data[i_var*ind.n_qpoint + ind.i_qpoint(i_row)];
        d = coef*d + w(i_row, i_var);
      }
    }
  }

  // reads from the face values associated with a given row
  static Bound read_bound(const double* face, Row_index ind)
  {
    Bound b;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        b(is_positive, i_var) = face[((ind.i_dim*2 + is_positive)*(ind.n_dim + 2) + i_var)*ind.n_fqpoint + ind.i_face_qpoint()];
      }
    }
    return b;
  }

  // writes the values in `b` to the face data associated with the specified row
  static void write_bound(Bound b, double* face, Row_index ind)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        face[((ind.i_dim*2 + is_positive)*(ind.n_dim + 2) + i_var)*ind.n_fqpoint + ind.i_face_qpoint()] = b(is_positive, i_var);
      }
    }
  }

  // reads from the face values associated with a given row
  static Bound read_bound(const std::array<double*, 6> faces, Row_index ind)
  {
    Bound b;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        b(is_positive, i_var) = faces[ind.i_dim*2 + is_positive][i_var*ind.n_fqpoint + ind.i_face_qpoint()];
      }
    }
    return b;
  }

  // writes the values in `b` to the face data associated with the specified row
  static void write_bound(Bound b, std::array<double*, 6> faces, Row_index ind)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        faces[ind.i_dim*2 + is_positive][i_var*ind.n_fqpoint + ind.i_face_qpoint()] = b(is_positive, i_var);
      }
    }
  }
};

}
#endif
