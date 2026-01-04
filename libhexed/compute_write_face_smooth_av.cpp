#include <hexed/kernel_utils.hpp>
#include <hexed/compute_write_face_smooth_av.hpp>

namespace hexed {

void compute_write_face_smooth_av(Kernel_mesh mesh) {
  (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis,
                                                                     mesh.n_var, 1., 1.))(mesh.elems);
}

}
