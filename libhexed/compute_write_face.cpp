#include <hexed/kernel_utils.hpp>
#include <hexed/compute_write_face.hpp>
#include <hexed/Navier_stokes.hpp>

namespace hexed {

typedef Navier_stokes<false, laminar> euler;
typedef Navier_stokes<false, k_omega> k_omega_euler;

void compute_write_face(Kernel_mesh mesh) {
  #define COMPUTE(pde_class) (*kernel_factory<Spatial<pde_class::Pde, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var))(mesh.elems);
  if (mesh.turb_model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

}
