#include <hexed/kernel_utils.hpp>
#include <hexed/compute_prolong_advection.hpp>

namespace hexed {

void compute_prolong_advection(Kernel_mesh mesh) {
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>>
    face_refs(mesh.face_refinements);
  (*kernel_factory<Spatial<pde::Advection, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis,
                                                                    mesh.mask_level, mesh.n_var, false, false))
                                                                   (face_refs);
}

}
