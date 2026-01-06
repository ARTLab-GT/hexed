#ifndef HEXED_COMPUTE_FACE_PERMUTATION_HPP_
#define HEXED_COMPUTE_FACE_PERMUTATION_HPP_

#include "Face_permutation.hpp"
#include "Transport_model.hpp"

namespace hexed {

std::unique_ptr<Face_permutation_dynamic> compute_face_permutation(int n_dim, int row_size, Connection_direction,
                                                                   double* data, Turbulence_model model);

}
#endif
