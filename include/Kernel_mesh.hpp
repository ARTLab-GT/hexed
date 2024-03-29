#ifndef HEXED_KERNEL_MESH_HPP_
#define HEXED_KERNEL_MESH_HPP_

#include "Basis.hpp"
#include "Sequence.hpp"
#include "Stopwatch_tree.hpp"
#include "Kernel_connection.hpp"
#include "Kernel_element.hpp"
#include "Refined_face.hpp"

namespace hexed
{

struct Kernel_mesh
{
  int n_dim;
  int row_size;
  const Basis& basis;
  Sequence<Kernel_connection&>& car_cons;
  Sequence<Kernel_connection&>& def_cons;
  Sequence<Kernel_element&>& car_elems;
  Sequence<Kernel_element&>& def_elems;
  Sequence<Kernel_element&>& elems;
  Sequence<Refined_face&>& ref_faces;
};

}
#endif
