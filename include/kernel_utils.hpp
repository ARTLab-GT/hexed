#ifndef HEXED_KERNEL_UTILS_HPP_
#define HEXED_KERNEL_UTILS_HPP_

#include <kernels.hpp>
#include <pde.hpp>
#include <Spatial.hpp>

#define COMPUTE_CONVECTION(Pde_templ, ...) { \
  Vector_view<Hard_kernel_connection> car_cons(mesh.car_connections); \
  Vector_view<Hard_kernel_connection> def_cons(mesh.def_connections); \
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs(mesh.face_refinements); \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.implicit_opts, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.implicit_opts, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
}

#define COMPUTE_DIFFUSION(Pde_templ, ...) { \
  Vector_view<Hard_kernel_connection> car_cons(mesh.car_connections); \
  Vector_view<Hard_kernel_connection> def_cons(mesh.def_connections); \
  Vector_view<std::vector<Kernel_face_refinement>&, std::vector<Kernel_face_refinement>> face_refs(mesh.face_refinements); \
  (*kernel_factory<Spatial<Pde_templ, false>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(car_cons, opts.sw_car, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor>(mesh.n_dim, mesh.row_size, opts.i_stage, mesh.mask_level, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(def_cons, opts.sw_def, "neighbor"); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, false, true))(face_refs, opts.sw_pr); \
  (*kernel_factory<Spatial<Pde_templ, false>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.implicit_opts, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "local"); \
  (*kernel_factory<Spatial<Pde_templ,  true>::Local>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, opts.implicit_opts, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "local"); \
  if (!opts.i_stage) { \
    (*kernel_factory<Spatial<Pde_templ, false>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, true, true))(face_refs, opts.sw_pr); \
    flux_bc(); \
    (*kernel_factory<Spatial<Pde_templ, false>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size, mesh.mask_level))(car_cons, opts.sw_car, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Neighbor_reconcile>(mesh.n_dim, mesh.row_size, mesh.mask_level))(def_cons, opts.sw_def, "neighbor"); \
    (*kernel_factory<Spatial<Pde_templ, false>::Restrict_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var, true, true))(face_refs, opts.sw_pr); \
    (*kernel_factory<Spatial<Pde_templ, false>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "reconcile LDG flux"); \
    (*kernel_factory<Spatial<Pde_templ,  true>::Reconcile_ldg_flux>(mesh.n_dim, mesh.row_size, mesh.basis, opts.dt, opts.i_stage, opts.compute_residual, opts.use_filter, opts.mask, opts.conv_substep, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "reconcile LDG flux"); \
  } \
  (*kernel_factory<Spatial<Pde_templ,  true>::Prolong_refined>(mesh.n_dim, mesh.row_size, mesh.basis, mesh.mask_level, mesh.n_var))(face_refs, opts.sw_pr); \
}

#define COMPUTE_MAX_DT(Pde_templ, ...) { \
  return std::min((*kernel_factory<Spatial<Pde_templ, false>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.car_elems, opts.sw_car, "compute time step"), \
                  (*kernel_factory<Spatial<Pde_templ,  true>::Max_dt>(mesh.n_dim, mesh.row_size, mesh.basis, local_time, opts.use_filter, msc, msd, calc_ts_ratio, mesh.n_var __VA_OPT__(,) __VA_ARGS__))(mesh.def_elems, opts.sw_def, "compute time step")); \
}

#endif
