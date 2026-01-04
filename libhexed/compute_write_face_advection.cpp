#include <hexed/kernel_utils.hpp>
#include <hexed/compute_write_face_advection.hpp>

namespace hexed {

void compute_write_face_advection(Kernel_mesh mesh, int offset) {
  (*kernel_factory<Spatial<pde::Advection, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var, 1.,
                                                               offset))(mesh.elems);
}

}
