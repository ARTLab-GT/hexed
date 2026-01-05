#include <hexed/kernel_utils.hpp>
#include <hexed/compute_prolong.hpp>
#include <hexed/Navier_stokes.hpp>

namespace hexed {

typedef Navier_stokes<false, laminar> euler;
typedef Navier_stokes<false, k_omega> k_omega_euler;
void compute_prolong(Kernel_mesh mesh, bool scale, bool offset) {
  #define COMPUTE(pde_class) { \
    Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs(mesh.face_refinements); \
    (*kernel_factory<Spatial<pde_class::Pde, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, scale, offset))(face_refs); \
  }
  if (mesh.turb_model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

}
