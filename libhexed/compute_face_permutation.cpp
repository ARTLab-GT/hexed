#include <hexed/kernel_utils.hpp>
#include <hexed/compute_face_permutation.hpp>

namespace hexed {

typedef pde::Navier_stokes<false, laminar> euler;
typedef pde::Navier_stokes<false, k_omega> k_omega_euler;

std::unique_ptr<Face_permutation_dynamic> compute_face_permutation(int n_dim, int row_size, Connection_direction dir,
                                                                   double* data, Turbulence_model model) {
  #define COMPUTE(pde_class) return kernel_factory<Spatial<pde_class::Pde, true>::Face_permutation>(n_dim, row_size, dir, data);
  if (model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

}
