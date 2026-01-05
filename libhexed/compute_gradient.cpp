#include <hexed/kernel_utils.hpp>
#include <hexed/compute_gradient.hpp>
#include <hexed/Gradient.hpp>

namespace hexed {

void compute_gradient(Kernel_mesh mesh, Kernel_options opts, std::function<void()> state_bc,
                      std::function<void()> flux_bc, int read_offset, int write_offset) {
  const int nq = math::pow(mesh.row_size, mesh.n_dim);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < mesh.elems.size(); ++i_elem) {
    double* tss = mesh.elems[i_elem].time_step_scale();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) tss[i_qpoint] = 1.;
  }
  (*kernel_factory<Spatial<Gradient, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var,
                                                         read_offset, write_offset))(mesh.elems);
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs
    (mesh.face_refinements);
  (*kernel_factory<Spatial<Gradient,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size,
                                                              mesh.basis, mesh.mask_level, mesh.n_var))
                                                             (face_refs, opts.sw_pr);
  state_bc();
  opts.dt = 1.;
  COMPUTE_DIFFUSION(Gradient, read_offset, write_offset)
}

}
