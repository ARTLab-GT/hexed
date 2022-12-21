#ifndef HEXED_ROW_RW_HPP_
#define HEXED_ROW_RW_HPP_

#include "Row_index.hpp"

namespace hexed
{

template <int n_var, int row_size>
class Row_rw
{
  public:
  typedef Mat<row_size, n_var> Row;
  typedef Mat<2, n_var> Bound;

  Row_rw() = delete;

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

  static void write_row(Row w, double* data, Row_index ind, double coef)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_row = 0; i_row < row_size; ++i_row) {
        double& d = data[i_var*ind.n_qpoint + ind.i_qpoint(i_row)];
        d = coef*d + w(i_row, i_var);
      }
    }
  }

  static Bound read_bound(const double* face, Row_index ind)
  {
    Bound b;
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        b(is_positive, i_var) = face[((ind.i_dim*2 + is_positive)*n_var + i_var)*ind.n_fqpoint + ind.i_face_qpoint()];
      }
    }
    return b;
  }

  static void write_bound(Bound b, double* face, Row_index ind)
  {
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int is_positive : {0, 1}) {
        face[((ind.i_dim*2 + is_positive)*n_var + i_var)*ind.n_fqpoint + ind.i_face_qpoint()] = b(is_positive, i_var);
      }
    }
  }
};

}
#endif
