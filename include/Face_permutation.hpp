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

};
#endif
