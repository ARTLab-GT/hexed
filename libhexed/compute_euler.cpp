#include <hexed/kernel_utils.hpp>
#include <hexed/compute_euler.hpp>
#include <hexed/Navier_stokes.hpp>
#include <hexed/Transport_model.hpp>

namespace hexed {

typedef Navier_stokes<false, laminar> euler;
typedef Navier_stokes<false, k_omega> k_omega_euler;

void compute_euler(Kernel_mesh mesh, Kernel_options opts) {
  if (mesh.turb_model == k_omega) COMPUTE_CONVECTION(k_omega_euler::Pde)
  else COMPUTE_CONVECTION(euler::Pde)
}

}
