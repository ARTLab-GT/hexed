#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed {

#define COMPUTE_CONVECTION(Pde_templ, ...) { \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(mesh.ref_faces, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.time_scheme, opts.be_dt, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(mesh.ref_faces, opts.sw_pr); \
}

typedef pde::Navier_stokes<false, laminar> euler;
typedef pde::Navier_stokes<false, k_omega> k_omega_euler;
void compute_euler(Kernel_mesh mesh, Kernel_options opts) {
  if (mesh.turb_model == k_omega) COMPUTE_CONVECTION(k_omega_euler::Pde)
  else COMPUTE_CONVECTION(euler::Pde)
}
void compute_advection(Kernel_mesh mesh, Kernel_options opts, double advect_length) COMPUTE_CONVECTION(pde::Advection, advect_length)

#undef COMPUTE_CONVECTION

void compute_prolong(Kernel_mesh mesh, bool scale, bool offset) {
  #define COMPUTE(pde_class) (*kernel_factory<Spatial<pde_class::Pde, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, scale, offset))(mesh.ref_faces);
  if (mesh.turb_model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

void compute_restrict(Kernel_mesh mesh, bool scale, bool offset) {
  #define COMPUTE(pde_class) (*kernel_factory<Spatial<pde_class::Pde, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, scale, offset))(mesh.ref_faces);
  if (mesh.turb_model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

void compute_prolong_advection(Kernel_mesh mesh) {
  (*kernel_factory<Spatial<pde::Advection, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, false, false))(mesh.ref_faces);
}

std::unique_ptr<Face_permutation_dynamic> face_permutation(int n_dim, int row_size, Connection_direction dir, double* data, Turbulence_model model) {
  #define COMPUTE(pde_class) return kernel_factory<Spatial<pde_class::Pde, true>::Face_permutation>(n_dim, row_size, dir, data);
  if (model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

void compute_write_face(Kernel_mesh mesh) {
  #define COMPUTE(pde_class) (*kernel_factory<Spatial<pde_class::Pde, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var))(mesh.elems);
  if (mesh.turb_model == k_omega) COMPUTE(k_omega_euler)
  else COMPUTE(euler)
  #undef COMPUTE
}

void compute_write_face_advection(Kernel_mesh mesh) {
  (*kernel_factory<Spatial<pde::Advection, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var, 1.))(mesh.elems);
}

void compute_write_face_smooth_av(Kernel_mesh mesh) {
  (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Write_face>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.n_var, 1., 1.))(mesh.elems);
}

}
