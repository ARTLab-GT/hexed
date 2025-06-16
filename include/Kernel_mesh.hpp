#ifndef HEXED_KERNEL_MESH_HPP_
#define HEXED_KERNEL_MESH_HPP_

#include "Basis.hpp"
#include "Sequence.hpp"
#include "Stopwatch_tree.hpp"
#include "Kernel_connection.hpp"
#include "Kernel_element.hpp"
#include "Refined_face.hpp"
#include "Transport_model.hpp"

namespace hexed {

struct Kernel_mesh {
  const int n_dim;
  const int row_size;
  const int n_var;
  const int mask_level;
  const Basis& basis;
  Turbulence_model turb_model; //! \note this is relevant because it affects the number/treatment of state variables
  Sequence<Kernel_element&>& car_elems;
  Sequence<Kernel_element&>& def_elems;
  Sequence<Kernel_element&>& elems;
  std::vector<Hard_kernel_connection> car_connections;
  std::vector<Hard_kernel_connection> def_connections;
  std::vector<std::vector<Kernel_face_refinement>> face_refinements;
};

}
#endif
