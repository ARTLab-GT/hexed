#include <hexed/kernel_utils.hpp>
#include <hexed/max_dt_eikonal.hpp>
#include <hexed/Eikonal.hpp>

namespace hexed {

double max_dt_eikonal(Kernel_mesh mesh, Kernel_options opts, double msc, double msd, double smoothing,
                      double grad_smoothing, double base_diff) {
  bool calc_ts_ratio = false;
  bool local_time = true;
  COMPUTE_MAX_DT(Eikonal, smoothing, grad_smoothing, base_diff)
  (*kernel_factory<Spatial<Eikonal, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var,
                                                        smoothing, grad_smoothing, base_diff))(mesh.elems);
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs
    (mesh.face_refinements);
  (*kernel_factory<Spatial<Eikonal,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size,
                                                             mesh.basis, mesh.mask_level, mesh.n_var))
    (face_refs, opts.sw_pr);
}

}
