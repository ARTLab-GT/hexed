#ifndef HEXED_FACE_PERMUTATION_HPP_
#define HEXED_FACE_PERMUTATION_HPP_

#include <Eigen/Dense>
#include "kernel_factory.hpp"
#include "Kernel_connection.hpp"
#include "math.hpp"

namespace hexed
{

/*
 * Permutes quadrature points on a face to match another face.
 * For deformed element connections, connecting across different dimensions can cause
 * the ordering of quadrature points on the respective faces to mismatch.
 * This class reorders the data corresponding to face 1 to match face 0.
 */
class Face_permutation_dynamic
{
  public:
  virtual ~Face_permutation_dynamic() = default;
  // reorder face 1 data to match face 0
  virtual void match_faces() = 0;
  // restore face 1 data to its original order
  virtual void restore() = 0;
};

template<int n_dim, int row_size>
class Face_permutation : public Face_permutation_dynamic
{
  static constexpr int n_qpoint = math::pow(row_size, n_dim - 1);
  static constexpr int n_var = n_dim + 2;
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

};
#endif
