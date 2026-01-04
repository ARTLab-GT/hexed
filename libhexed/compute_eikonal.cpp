#include <hexed/kernel_utils.hpp>
#include <hexed/compute_eikonal.hpp>

namespace hexed {

void compute_eikonal(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, std::function<void()> state_bc,
                     std::function<void()> flux_bc, double smoothing, double grad_smoothing, double base_diff) {
  const int nq = math::pow(mesh.row_size, mesh.n_dim);
  double dt = opts.dt;
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs
    (mesh.face_refinements);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < mesh.elems.size(); ++i_elem) {
    double* tss = mesh.elems[i_elem].time_step_scale();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) tss[i_qpoint] = 1.;
  }
  const int read_offset = pde::laplacian_av_offset(mesh.n_var);
  const int write_offset = pde::residual_cache_offset(mesh.n_var, mesh.row_size) + 1 + mesh.n_dim;
  (*kernel_factory<Spatial<pde::Laplace, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var,
                                                             read_offset, write_offset))(mesh.elems);
  (*kernel_factory<Spatial<pde::Laplace,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size,
                                                                  mesh.basis, mesh.mask_level, mesh.n_var))
                                                                 (face_refs, opts.sw_pr);
  state_bc();
  opts.dt = 1.;
  COMPUTE_DIFFUSION(pde::Laplace, read_offset, write_offset)
  opts.dt = dt;
  max_dt_eikonal(mesh, opts, msc, msd, smoothing, grad_smoothing, base_diff);
  (*kernel_factory<Spatial<pde::Eikonal, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var,
                                                             smoothing, grad_smoothing, base_diff))(mesh.elems);
  (*kernel_factory<Spatial<pde::Eikonal,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size,
                                                                  mesh.basis, mesh.mask_level, mesh.n_var))
                                                                 (face_refs, opts.sw_pr);
  state_bc();
  COMPUTE_DIFFUSION(pde::Eikonal, smoothing, grad_smoothing, base_diff)
}

}
