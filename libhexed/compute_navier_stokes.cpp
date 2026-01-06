#include <hexed/kernel_utils.hpp>
#include <hexed/compute_navier_stokes.hpp>
#include <hexed/Navier_stokes.hpp>

namespace hexed {

typedef Navier_stokes<true, laminar> ns;
typedef Navier_stokes<true, k_omega> rans;
void compute_navier_stokes(Kernel_mesh mesh, Kernel_options opts, std::function<void()> flux_bc,
                           Transport_model visc, Transport_model therm_cond) {
  if (mesh.turb_model == k_omega) COMPUTE_DIFFUSION(rans::Pde, visc, therm_cond)
  else COMPUTE_DIFFUSION(ns::Pde, visc, therm_cond)
}

}
