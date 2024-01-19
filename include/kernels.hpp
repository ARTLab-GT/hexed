#ifndef HEXED_KERNELS_HPP_
#define HEXED_KERNELS_HPP_

#include "Basis.hpp"
#include "Sequence.hpp"
#include "Kernel_connection.hpp"
#include "Kernel_element.hpp"
#include "Refined_face.hpp"
#include "Stopwatch_tree.hpp"

namespace hexed
{

struct Kernel_args
{
  int n_dim;
  int row_size;
  const Basis& basis;
  Stopwatch_tree& sw_car;
  Stopwatch_tree& sw_def;
  Stopwatch_tree& sw_pr;
  Sequence<Kernel_connection&>& car_cons;
  Sequence<Kernel_connection&>& def_cons;
  Sequence<Kernel_element&>& car_elems;
  Sequence<Kernel_element&>& def_elems;
  Sequence<Refined_face&>& ref_faces;
  double dt;
  int i_stage;
  bool compute_residual = false;
  bool use_filter = false;
};

void compute_euler(Kernel_args&);
void compute_advection(Kernel_args&);

}
#endif
