#include <fstream>
#include <iostream>
#include <H5Cpp.h>

#include <hexed/config.hpp>
#include <hexed/Solver.hpp>
#include <hexed/Tecplot_file.hpp>
#include <hexed/Vis_data.hpp>
#include <hexed/Xdmf_wrapper.hpp>
#include <hexed/iterative.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include <hexed/Face_permutation.hpp>
#include <hexed/Row_index.hpp>
#include <hexed/stabilizing_art_visc.hpp>
#include <hexed/Array.hpp>
#include <hexed/vis_variables.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

const double dirk2_gamma = 1 - std::sqrt(.5);

Kernel_mesh& Solver::_kernel_mesh() {
  return _preti_masks[0]->kernel_mesh;
}

void Solver::_put_cache() {
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    double* cache = elems[i_elem].residual_cache();
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) cache[i_dof] = state[i_dof];
  }
}

void Solver::_get_cache() {
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    double* cache = elems[i_elem].residual_cache();
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = cache[i_dof];
  }
}

double max_fun(double x, double y) {return std::max(x, y);}
double min_fun(double x, double y) {return std::min(x, y);}

void Solver::share_vertex_data(std::function<double&(Element&, int i_vertex)> access_fun, bool minmax) {
  share_vertex_data(access_fun, access_fun, minmax);
}

void Solver::share_vertex_data(std::function<double(Element&, int i_vertex)> get,
                               std::function<double&(Element&, int i_vertex)> set,
                               bool minmax) {
  int nv = params.n_vertices();
  auto verts = acc_mesh->shape_vertices();
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (Int i_vert = 0; i_vert < (Int)verts.size(); ++i_vert) {
    next::Vertex::Shared_value(verts[i_vert]).set(minmax ? -huge : huge);
  }
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    auto& elem = elements[i_elem];
    auto& shape = elem.shape();
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      next::Vertex::Shared_value shared(shape.vertex(i_vert));
      shared.set(get(elem, i_vert), minmax);
    }
  }
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    auto& elem = elements[i_elem];
    auto& shape = elem.shape();
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      next::Vertex::Shared_value shared(shape.vertex(i_vert));
      set(elem, i_vert) = shared.get();
    }
  }
}

void Solver::apply_state_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto bc_cons {_preti_masks[_preti_level]->bound_cons};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con]->boundary_condition();
    acc_mesh->boundary_condition(bc_sn).apply_state(*bc_cons[i_con]);
    Array<double> ghost = bc_cons[i_con]->ghost().flow_state()(0);
    Array<double> inside = bc_cons[i_con]->inside().flow_state()(0);
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < params.n_face_qpoint(); ++i_qpoint) {
        if (!(std::abs(ghost(i_var)[i_qpoint]) < 1e10)) ghost(i_var)[i_qpoint] = inside(i_var)[i_qpoint];
      }
    }
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_flux_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto bc_cons {_preti_masks[_preti_level]->bound_cons};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    Boundary_connection* con = bc_cons[i_con];
    con->flux_cache() = con->inside().flow_state()(1);
    int bc_sn = con->boundary_condition();
    acc_mesh->boundary_condition(bc_sn).apply_flux(*con);
    Array<double> ghost = bc_cons[i_con]->ghost().flow_state()(1);
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < params.n_face_qpoint(); ++i_qpoint) {
        if (!(std::abs(ghost(i_var)[i_qpoint]) < 1e10)) ghost(i_var)[i_qpoint] = 0.;
      }
    }
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_avc_diff_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].boundary_condition();
    acc_mesh->boundary_condition(bc_sn).apply_diffusion(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_avc_diff_flux_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].boundary_condition();
    acc_mesh->boundary_condition(bc_sn).flux_diffusion(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_fta_flux_bcs() {
  auto bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    bc_cons[i_con].ghost().flow_state()(1) = -bc_cons[i_con].inside().flow_state()(1);
  }
}

bool Solver::use_ldg() {
  return visc.is_viscous || therm_cond.is_viscous || use_art_visc;
}

double Solver::max_dt(double msc, double msd, double lim_thresh) {
  Kernel_options opts {
    stopwatch["cartesian"],
    stopwatch["deformed"],
    stopwatch["prolong/restrict"],
    0, 0, bool(_namespace->get<int>("use_filter")),
  };
  bool local_time = _time_scheme != explicit_unsteady;
  if (use_ldg()) {
    return max_dt_navier_stokes(_kernel_mesh(), opts, msc, msd, lim_thresh, local_time, visc, therm_cond);
  } else {
    return max_dt_euler(_kernel_mesh(), opts, msc, msd, local_time);
  }
}

void Solver::_init_face_state() {
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  auto bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].boundary_condition();
    acc_mesh->boundary_condition(bc_sn).init_cache(bc_cons[i_con]);
  }
  update_bound_conds();
}

Interpreter Solver::_interpreter() {
  Interpreter inter(std::vector<std::string>{});
  inter.variables = _namespace;
  return inter;
}

Solver::Solver(int n_dim, int row_size, double root_mesh_size, Time_scheme time_scheme,
               Transport_model viscosity_model, Transport_model thermal_conductivity_model,
               Turbulence_model turbulence_model,
               std::shared_ptr<Namespace> space)
: params{2 + n_extra_stage(time_scheme), n_dim + 2 + 2*(turbulence_model == k_omega), n_dim, row_size}
, acc_mesh{new Accessible_mesh(params, root_mesh_size, turbulence_model)}
, basis{row_size}
, stopwatch{"(element*update)"}
, use_art_visc{false}
, fix_admis{false}
, av_rs{row_size}
, visc{viscosity_model}
, therm_cond{thermal_conductivity_model}
, turb{turbulence_model}
, _namespace{space}
, _implicit{false}
, _preti_level{0}
, _time_scheme{time_scheme}
{
  _namespace->assign_default("max_safety", .7); // maximum allowed safety factor for time stepping
  _namespace->assign_default("max_time_step", huge); // maximum allowed time step
  _namespace->assign_default("fix_admis_max_safety", .2); // staility ratio for fixing thermodynamic admissibility.
  _namespace->assign_default("av_diff_ratio", .3); // ratio of diffusion time to advection width
  // final scaling parameter applied to artificial viscosity coefficient
  _namespace->assign_default("av_visc_mult", 1e2);
  // maximum artificial viscosity coefficient before scaling (i.e. nondimensional)
  _namespace->assign_default("av_unscaled_max", 5.);
  _namespace->assign_default("av_advect_max_safety", .7); // stability ratio for advection
  _namespace->assign_default("av_diff_max_safety", .7); // stability ratio for diffusion
  _namespace->assign_default("vis_art_visc_vars", 0);
  _namespace->assign_default("n_cheby_bl", 1);
  _namespace->assign_default("n_cheby_flow", 1);
  _namespace->assign_default("max_conv_sub_iters", 1);
  _namespace->assign_default("n_cheby_av", 1);
  _namespace->assign_default("cheby_safety", .9); // safety factor to apply to Chebyshev-acceleration
  _namespace->assign_default("bl_multirate", 0);
  // number of advection iterations to run each time `update_art_visc_smoothness` is called
  _namespace->assign_default("av_advect_iters", 1);
  // number of diffusion iterations to run each time `update_art_visc_smoothness` is called
  _namespace->assign_default("av_diff_iters", 1);
  _namespace->assign_default("flow_iters", 1);
  _namespace->assign_default("bl_iters", 1);
  _namespace->assign_default("fix_iters", 0);
  _namespace->assign_default("use_filter", 0); // whether to use modal filter acceleration
  _namespace->assign_default("elementwise_art_visc", 0);
  _namespace->assign_default("elementwise_art_visc_diff_ratio", 5.);
  _namespace->assign_default<std::string>("working_dir", ".");
  _namespace->assign_default("iteration", 0);
  _namespace->assign_default("pseudotime_iteration", 0);
  _namespace->assign_default("flow_time", 0.);
  _namespace->assign_default("geom_length", 0.);
  if (!is_implicit(_time_scheme)) _namespace->assign_default("time_step", 0.);
  _namespace->assign("time_stage", 0);
  _namespace->assign("n_time_stages", n_total_stage(_time_scheme));
  _namespace->assign_default("art_visc_residual", 0.);
  status.set_time();
  // setup categories for performance reporting
  std::string unit = "(element*(time integration stage))";
  stopwatch.emplace("prolong/restrict", unit);
  stopwatch.emplace("fix admis.", "(element*(fix admis. iter))");
  stopwatch["fix admis."].emplace("check admis.", "(element*update)");
  stopwatch.emplace("set art visc", stopwatch.work_unit_name);
  stopwatch["set art visc"].emplace("initialize", stopwatch.work_unit_name);
  stopwatch["set art visc"].emplace("advection", stopwatch.work_unit_name);
  stopwatch["set art visc"]["advection"].emplace("update", unit);
  stopwatch["set art visc"]["advection"].emplace("setup", stopwatch.work_unit_name);
  stopwatch["set art visc"]["advection"].emplace("BCs", unit);
  stopwatch["set art visc"].emplace("diffusion", stopwatch.work_unit_name);
  for (std::string type : {"cartesian", "deformed"}) {
    for (auto* sw : {&stopwatch, &stopwatch["set art visc"]["advection"],
                     &stopwatch["set art visc"]["diffusion"], &stopwatch["fix admis."]}) {
      sw->emplace(type, stopwatch.work_unit_name);
      (*sw)[type].emplace("compute time step", stopwatch.work_unit_name);
      (*sw)[type].emplace("neighbor", "(connection*(time integration stage))");
      (*sw)[type].emplace("local", unit);
    }
    for (auto* sw : {&stopwatch, &stopwatch["fix admis."], &stopwatch["set art visc"]["diffusion"]}) {
      (*sw)[type].emplace("reconcile LDG flux", unit);
    }
  }
  stopwatch.emplace("boundary conditions", "(boundary connection)*(time integration stage)");
  stopwatch.emplace("visualization", "file");
  stopwatch["visualization"].emplace("field", "element");
  stopwatch["visualization"].emplace("field wireframe", "element");
  stopwatch["visualization"].emplace("surface", "surface face");
  stopwatch["visualization"].emplace("surface wireframe", "surface face");
  stopwatch["visualization"].emplace("contour", "element");
  stopwatch.emplace("integrals", "integral");
  stopwatch["integrals"].emplace("field", "integral");
  stopwatch["integrals"].emplace("surface", "integral");
  // initialize advection state to 1
  auto& elements = acc_mesh->elements();
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* adv = elements[i_elem].advection_state();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      for (int i_adv = 0; i_adv < rs; ++i_adv) adv[i_adv*nq + i_qpoint] = 1.;
    }
  }
}

Namespace& Solver::nspace() {return *_namespace;}

Mesh& Solver::mesh() {return *acc_mesh;}
Storage_params Solver::storage_params() {return params;}
const Stopwatch_tree& Solver::stopwatch_tree() {return stopwatch;}

void Solver::read_mesh(std::string file_name, std::vector<Flow_bc*> extremal_bcs,
                       Surface_geom* geom, Flow_bc* surface_bc) {
  acc_mesh.reset(new Accessible_mesh(file_name, extremal_bcs, turb, geom, surface_bc));
  HEXED_ASSERT(acc_mesh->storage_params().n_stage == params.n_stage,
               "attempt to read a mesh file with a different `n_stage`");
  HEXED_ASSERT(acc_mesh->storage_params().n_var == params.n_var,
               "attempt to read a mesh file with a different `n_var`");
  HEXED_ASSERT(acc_mesh->storage_params().n_dim == params.n_dim,
               "attempt to read a mesh file with a different `n_dim`");
  HEXED_ASSERT(acc_mesh->storage_params().row_size == params.row_size,
               "attempt to read a mesh file with a different `row_size`");
  HEXED_ASSERT(acc_mesh->storage_params().n_forcing == params.n_forcing,
               "attempt to read a mesh file with a different `n_forcing`");
  calc_jacobian();
}

void Solver::read_state(std::string file_name) {
  auto& elems = acc_mesh->elements();
  H5::H5File file(file_name + ".state.h5", H5F_ACC_RDONLY);
  hsize_t n_elem = elems.size();
  hsize_t n_var = params.n_var_numeric();
  hsize_t n_qpoint = params.n_qpoint();
  auto dset = file.openDataSet("state");
  for (hsize_t i_elem = 0; i_elem < n_elem; ++i_elem) {
    hsize_t elem_dims [3] {1, n_var, n_qpoint};
    H5::DataSpace mspace (3, elem_dims, nullptr);
    hsize_t offset [3] {i_elem, 0, 0};
    hsize_t stride [3] {1, 1, 1};
    hsize_t block [3] {1, 1, 1};
    auto dspace = dset.getSpace();
    dspace.selectHyperslab(H5S_SELECT_SET, elem_dims, offset, stride, block);
    dset.read(elems[i_elem].state(), dset.getDataType(), mspace, dspace);
  }
  _init_face_state();
}

void Solver::write_state(std::string file_name) {
  auto& elems = acc_mesh->elements();
  H5::H5File file(file_name + ".state.h5", H5F_ACC_TRUNC);
  hsize_t n_elem = elems.size();
  hsize_t n_var = params.n_var_numeric();
  hsize_t n_qpoint = params.n_qpoint();
  hsize_t dims [3] {n_elem, n_var, n_qpoint};
  H5::DataSpace dspace(3, dims);
  auto dset = file.createDataSet("state", H5::PredType::NATIVE_DOUBLE, dspace);
  for (hsize_t i_elem = 0; i_elem < n_elem; ++i_elem) {
    hsize_t elem_dims [3] {1, n_var, n_qpoint};
    H5::DataSpace mspace (3, elem_dims, nullptr);
    hsize_t offset [3] {i_elem, 0, 0};
    hsize_t stride [3] {1, 1, 1};
    hsize_t block [3] {1, 1, 1};
    dspace.selectHyperslab(H5S_SELECT_SET, elem_dims, offset, stride, block);
    dset.write(elems[i_elem].state(), H5::PredType::NATIVE_DOUBLE, mspace, dspace);
  }
}

void Solver::calc_jacobian() {
  acc_mesh->valid().assert_valid();
  _preti_masks = acc_mesh->preti_masks(basis, _namespace->get<int>("preti"));
  // compute element jacobians
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].set_jacobian(basis);
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
      //HEXED_ASSERT(elements[i_elem].jacobian_determinant(i_qpoint) > 0., "Nonpositive Jacobian")
      double det = elements[i_elem].jacobian_determinant(i_qpoint);
      if (!(det > 0. && std::isfinite(det))) {
        std::string message = "Nonpositive Jacobian (" + to_string(det) + "). Node positions:\n";
        for (int i_row = 0; i_row < params.row_size; ++i_row) {
          for (int j_row = 0; j_row < params.row_size; ++j_row) {
            for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
              message += to_string(elements[i_elem].position(basis)(i_dim)[i_row*params.row_size + j_row]) + ",";
            }
            message += "  ";
          }
          message += "\n";
        }
        elements[i_elem].shape().visualize("default", "bad_element");
        HEXED_THROW(message)
      }
    }
  }
  // do some extra work to make sure each face knows its normal vectors
  auto face_refs = acc_mesh->face_refinements();
  #pragma omp parallel for
  for (auto& vec : face_refs) {
    for (auto& ref : vec) if (ref.is_deformed()) {
      ref.coarse().flow_state()(1)(0, params.n_dim) = ref.coarse().normal();
    }
  }
  compute_prolong(_preti_masks[0]->kernel_mesh, 0, 1);
  #pragma omp parallel for
  for (auto& vec : face_refs) {
    for (auto& ref : vec) if (ref.is_deformed()) {
      for (int i_fine = 0; i_fine < 2; ++i_fine) {
        ref.fine()[i_fine]->normal() = ref.fine()[i_fine]->flow_state()(1)(0, params.n_dim);
      }
    }
  }
  // set position at boundary faces
  auto bc_cons = acc_mesh->boundary_connections();
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    HEXED_ASSERT(con.inside().element(), "connection has no element")
    Element& elem = *con.inside().element();
    con.position() = elem.face_position(basis)(con.inside().i_dim())(con.inside().sign());
    if (con.inside().is_deformed()) con.ghost().normal() = con.inside().normal();
  }
  share_vertex_data(&Element::vertex_time_step_scale, false);
  // check that all the normals agree on both faces of every connection
  for (Neighbor_connection& con : acc_mesh->neighbor_connections(1)) {
    auto dir = con.get_direction();
    Array<double> temp_storage = Array<double>::make_uniform({params.n_var, params.n_face_qpoint()}, 0.);
    temp_storage(0, params.n_dim) = con.face(1).normal();
    auto perm = face_permutation(params.n_dim, params.row_size, dir, temp_storage.data(), turb);
    perm->match_faces();
    Array<double> nrml0 = con.face(0).normal()*(con.face(0).nominal_area()*math::sign(!dir.flip_normal(0)));
    Array<double> nrml1 = temp_storage(0, params.n_dim)*(con.face(1).nominal_area()*math::sign(!dir.flip_normal(1)));
    HEXED_ASSERT((nrml0 - nrml1).norm() < 1e-3*con.face(0).nominal_area(),
                 str_cat("normal mismatch: ", dir, "\n", nrml0, nrml1))
  }
  _preti_masks[0]->desired_iters = 1;
  _preti_masks[0]->repeat = true;
  _init_face_state();
}

void Solver::initialize(std::string(expr)) {
  acc_mesh->assert_valid();
  std::vector<std::string> state_vars;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) state_vars.push_back("momentum" + std::to_string(i_dim));
  state_vars.push_back("density");
  state_vars.push_back("energy");
  if (turb == k_omega) {
    state_vars.push_back("turbulent_kinetic_energy");
    state_vars.push_back("turbulent_dissipation_bassi");
  }
  auto inter = _interpreter();
  int n_var = state_vars.size();
  int nq = params.n_qpoint();
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    auto& elem = elements[i_elem];
    auto sub = inter.make_sub();
    vis_variables::element(*sub.variables, elem);
    vis_variables::position(*sub.variables, elem, basis);
    sub.exec(expr);
    Array<double> state({n_var, nq}, elem.state());
    for (int i_var = 0; i_var < n_var; ++i_var) {
      sub.variables->assign_array(state(i_var), state_vars[i_var]);
    }
    for (int i_adv = 0; i_adv < params.n_advection(params.row_size); ++i_adv) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        elem.advection_state()[i_adv*nq + i_qpoint] = 1.;
      }
    }
    Array<double>({n_var, nq}, elem.residual_cache()) = 0;
  }
  if (is_implicit(_time_scheme)) _init_stage_storage(0);
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    for (int i_face = 0; i_face < 2*params.n_dim; ++i_face) {
      elems[i_elem].face(i_face).flow_state()(1) = 0; // initialize face flux to 0
    }
  }
  _init_face_state();
}

bool Solver::using_art_visc() {
  return use_art_visc;
}

void Solver::set_art_visc_off() {
  use_art_visc = false;
}

void Solver::set_art_visc_constant(double value) {
  use_art_visc = true;
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* av = elements[i_elem].bulk_av_coef();
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
     av[i_qpoint] = value;
    }
  }
}

void Solver::diffuse_art_visc(double diff_time) {
  // evaluate CFL condition
  double diff_safety = _namespace->get<double>("av_diff_max_safety");
  double n_cheby = _namespace->get<double>("n_cheby_av");
  double cheby_safety = _namespace->get<double>("cheby_safety");
  Kernel_options opts {
    stopwatch["set art visc"]["diffusion"]["cartesian"],
    stopwatch["set art visc"]["diffusion"]["deformed"],
    stopwatch["prolong/restrict"],
    0.,
    0,
    false,
    false,
  };
  max_dt_smooth_av(_kernel_mesh(), opts, 1., diff_safety, true);
  // initialize residual to zero (will compute RMS over all real time steps)
  compute_write_face_smooth_av(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  // perform pseudotime iteration
  for (int i_iter = 0; i_iter < _namespace->get<int>("av_diff_iters"); ++i_iter) {
    HEXED_ASSERT(_preti_masks.size(), "meshing mask list is empty");
    int n_preti = (_namespace->get<int>("bl_multirate") && !i_iter) ? _preti_masks.size() : 1;
    for (int i_preti = 0; i_preti < n_preti; ++i_preti) {
      _preti_level = i_preti;
      int n_bl = i_preti ? _namespace->get<int>("bl_iters") : 1;
      Kernel_mesh& km = _preti_masks[_preti_level]->kernel_mesh;
      opts.mask = i_preti;
      for (int i_bl = 0; i_bl < n_bl; ++i_bl) {
        for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby) {
          double s = math::chebyshev_step(n_cheby, i_cheby, cheby_safety);
          apply_avc_diff_bcs();
          opts.dt = s;
          compute_smooth_av(km, opts, [this](){apply_avc_diff_flux_bcs();}, diff_time, s);
        }
      }
    }
  }
}

void Solver::update_art_visc_smoothness(double advect_length) {
  stopwatch.stopwatch.start();
  stopwatch["set art visc"].stopwatch.start();
  use_art_visc = true;
  const int nq = params.n_qpoint();
  const int nd = params.n_dim;
  const int rs = params.row_size;
  double heat_rat = 1.4;
  auto& elements = acc_mesh->elements();

  stopwatch["set art visc"]["initialize"].stopwatch.start();
  // set advection velocity
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].state();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double scale = sqrt(2*state[nd*nq + i_qpoint]*state[(nd + 1)*nq + i_qpoint]);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        state[i_dim*nq + i_qpoint] /= scale;
      }
    }
  }
  stopwatch["set art visc"]["initialize"].stopwatch.pause();
  stopwatch["set art visc"]["initialize"].work_units_completed += elements.size();
  // enforce CFL condition
  auto& sw_adv = stopwatch["set art visc"]["advection"];
  sw_adv.stopwatch.start();
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  double adv_safety = _namespace->get<double>("av_advect_max_safety");
  Kernel_options opts {
    sw_adv["cartesian"],
    sw_adv["deformed"],
    stopwatch["prolong/restrict"],
    1.,
    0,
    false,
    false,
  };
  max_dt_advection(_kernel_mesh(), opts, adv_safety, 1., true, advect_length);
  compute_write_face_advection(_kernel_mesh());
  compute_prolong_advection(_kernel_mesh());

  // begin estimation of high-order derivative in the style of the Cauchy-Kovalevskaya theorem
  // using a linear advection equation.

  // perform pseudotime iteration
  for (int iter = 0; iter < _namespace->get<int>("av_advect_iters"); ++iter) {
    HEXED_ASSERT(_preti_masks.size(), "meshing mask list is empty");
    int n_preti = (_namespace->get<int>("bl_multirate") && !iter) ? _preti_masks.size() : 1;
    for (int i_preti = 0; i_preti < n_preti; ++i_preti) {
      _preti_level = i_preti;
      int n_bl = i_preti ? _namespace->get<int>("bl_iters") : 1;
      Kernel_mesh& km = _preti_masks[_preti_level]->kernel_mesh;
      opts.mask = i_preti;
      for (int i_bl = 0; i_bl < n_bl; ++i_bl) {
        sw_adv["setup"].stopwatch.start();
        // evaluate advection operator
        sw_adv["setup"].stopwatch.pause();
        for (int i = 0; i < 2; ++i) {
          sw_adv["BCs"].stopwatch.start();
          auto bc_cons {acc_mesh->boundary_connections()};
          #pragma omp parallel for
          for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
            int bc_sn = bc_cons[i_con].boundary_condition();
            acc_mesh->boundary_condition(bc_sn).apply_advection(bc_cons[i_con]);
          }
          sw_adv["BCs"].stopwatch.pause();
          sw_adv["BCs"].work_units_completed += acc_mesh->elements().size();
          opts.i_stage = i;
          compute_advection(km, opts, advect_length);
        }
      }
    }
    sw_adv["cartesian"].work_units_completed += acc_mesh->cartesian().elements().size();
    sw_adv["deformed" ].work_units_completed += acc_mesh->deformed ().elements().size();
  }
  sw_adv["setup"].work_units_completed += elements.size();
  sw_adv["update"].work_units_completed += elements.size();
  stopwatch["set art visc"]["advection"].stopwatch.pause();
  stopwatch["set art visc"]["advection"].work_units_completed += elements.size();
  // compute projection onto Legendre polynomial
  Eigen::VectorXd weights = basis.node_weights();
  Eigen::VectorXd orth = basis.orthogonal(av_rs - 1);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* forcing = elements[i_elem].art_visc_forcing();
    double* adv = elements[i_elem].advection_state();
    double* state = elements[i_elem].state();
    double has_shock = false;
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double proj = 0;
      for (int i_proj = 0; i_proj < rs; ++i_proj) {
        proj += adv[i_proj*nq + i_qpoint]*weights(i_proj)*orth(i_proj);
      }
      double mach_suppression = 0;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        mach_suppression += state[i_dim*nq + i_qpoint]*state[i_dim*nq + i_qpoint];
      }
      mach_suppression /= heat_rat*(heat_rat - 1.);
      mach_suppression = mach_suppression*mach_suppression/(.3 + mach_suppression*mach_suppression);
      double f = proj*proj*2*state[(nd + 1)*nq + i_qpoint]/state[nd*nq + i_qpoint]*mach_suppression;
      forcing[i_qpoint] = std::isfinite(f) ? std::max(0., std::min(f, 1e10*advect_length*advect_length)) : 0.;
      has_shock = has_shock || forcing[i_qpoint] > 5.*advect_length*advect_length;
    }
    elements[i_elem].has_shock = has_shock;
    elements[i_elem].spread_shock = false;
  } // Cauchy-Kovalevskaya-style derivative estimate complete!

  for (int spread_iter = 0; spread_iter < 0; ++spread_iter) {
    for (bool is_def : {0, 1}) {
      #pragma omp parallel for
      for (Neighbor_connection& con : acc_mesh->neighbor_connections(is_def)) {
        bool shock = false;
        bool has_elems = true;
        for (int i_side = 0; i_side < 2; ++i_side) {
          Element* elem = con.face(i_side).find_element();
          has_elems = has_elems && elem;
          if (elem) shock = shock || elem->has_shock;
        }
        if (!shock || !has_elems) continue;
        for (int i_side = 0; i_side < 2; ++i_side) {
          #pragma omp atomic write
          con.face(i_side).find_element()->spread_shock = true;
        }
      }
    }
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      auto& elem = elements[i_elem];
      elem.has_shock = elem.has_shock || elem.spread_shock;
    }
  }

  // begin root-smear-square operation
  int n_real = params.n_forcing - 1; // number of real time steps (as apposed to pseudotime steps)
  // compute size of real time step (as opposed to pseudotime)
  double diff_time = _namespace->get<double>("av_diff_ratio")*advect_length*advect_length/n_real;
  stopwatch["set art visc"]["diffusion"].stopwatch.start();
  diffuse_art_visc(diff_time);
  stopwatch["set art visc"]["diffusion"].stopwatch.pause();
  stopwatch["set art visc"]["diffusion"].work_units_completed += elements.size();

  // clean up
  double mult = _namespace->get<double>("av_visc_mult")*advect_length;
  double us_max = advect_length*_namespace->get<double>("av_unscaled_max")
                  *std::sqrt(2*_namespace->get<double>("freestream" + std::to_string(nd + 1))
                  /_namespace->get<double>("freestream" + std::to_string(nd)));
  double resid = 0;
  Mat<> qpoint_weights = math::pow_outer(basis.node_weights(), nd);
  #pragma omp parallel for reduction(+:resid)
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].state();
    double* av = elements[i_elem].bulk_av_coef();
    double* forcing = elements[i_elem].art_visc_forcing();
    double volume = math::pow(elements[i_elem].nominal_size(), nd);
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double f = mult*forcing[n_real*nq + i_qpoint];
      double new_av = us_max*f/(us_max + f);
      resid += math::pow(av[i_qpoint] - new_av, 2)*qpoint_weights(i_qpoint)*volume;
      av[i_qpoint] = new_av;
      // put the flow state back how we found it
      double scale = sqrt(2*state[nd*nq + i_qpoint]*state[(nd + 1)*nq + i_qpoint]);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        state[i_dim*nq + i_qpoint] *= scale;
      }
    }
  }
  _namespace->assign("art_visc_residual", std::sqrt(resid));
  // update the face state
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  stopwatch["set art visc"].stopwatch.pause();
  stopwatch["set art visc"].work_units_completed += elements.size();
  stopwatch.stopwatch.pause();
}

void Solver::_init_stage_storage(int stage) {
  int n_var = params.n_var;
  int nq = params.n_qpoint();
  double time_step = _namespace->get<double>("time_step");
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    Array<double> state({n_var, nq}, elem.state());
    int n_res_cache = 1 + elem.get_is_deformed() + n_extra_stage(_time_scheme);
    Array<double> res_cache({n_res_cache, n_var, nq}, elem.residual_cache());
    Array<double> tss({nq}, elem.time_step_scale());
    if (_time_scheme == backward_euler) {
      res_cache(n_res_cache - 1) = state/time_step;
    } else if (_time_scheme == crank_nicolson) {
      for (int i_var = 0; i_var < n_var; ++i_var) {
        res_cache(n_res_cache - 1)(i_var) = res_cache(0)(i_var)/tss + state(i_var)/(.5*time_step);
      }
    } else if (_time_scheme == dirk2) {
      if (stage) {
        for (int i_var = 0; i_var < n_var; ++i_var) {
          res_cache(n_res_cache - 1)(i_var) += res_cache(0)(i_var)/tss*(1 - dirk2_gamma)/dirk2_gamma;
        }
      } else {
        res_cache(n_res_cache - 1) = state/(dirk2_gamma*time_step);
      }
    }
  }
}

int Solver::next_time_stage() {
  HEXED_ASSERT(is_implicit(_time_scheme), "This function is only for implicit time integration.")
  int stage = (_namespace->get<int>("time_stage") + 1)%n_total_stage(_time_scheme);
  if (_time_scheme != backward_euler) compute_residual();
  _init_stage_storage(stage);
  return stage;
}

void Solver::update_art_visc_elwise(double width, bool pde_based) {
  use_art_visc = true;
  Mass mass;
  set_uncertainty(Normalized_nonsmooth(mass));
  auto& elems = acc_mesh->elements();
  double scale = width/(basis.row_size - 1)*(_namespace->get<double>("freestream_speed") + _namespace->get<double>("freestream_sound_speed"));
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double& u = elems[i_elem].uncertainty;
    u = 2*std::log(u)/std::log(10);
    double ramp_center = -4 - 4.25*std::log(basis.row_size - 1)/std::log(10);
    double half_width = 0.5;
    if (!(u > ramp_center - half_width)) u = 0;
    else if (u >= ramp_center + half_width) u = 1;
    else u = .5*(1 + std::sin(constants::pi*(u - ramp_center)/2/half_width));
    u *= scale;
  }
  if (pde_based) {
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      double elem_av = elems[i_elem].uncertainty;
      double* av = elems[i_elem].laplacian_av_coef();
      double* forcing = elems[i_elem].art_visc_forcing();
      for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
        forcing[i_qpoint] = elem_av;
        forcing[params.n_qpoint() + i_qpoint] = av[i_qpoint];
      }
    }
    diffuse_art_visc(_namespace->get<double>("elementwise_art_visc_diff_ratio")*width*width);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      double* av = elems[i_elem].laplacian_av_coef();
      double* forcing = elems[i_elem].art_visc_forcing();
      for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) av[i_qpoint] = forcing[params.n_qpoint() + i_qpoint];
    }
    compute_write_face(_kernel_mesh());
    compute_prolong(_kernel_mesh());
  } else {
    share_vertex_data([](Element& elem, int){return elem.uncertainty;},
                      [](Element& elem, int i_vert)->double&{return elem.vertex_elwise_av(i_vert);},
                      true);
    Mat<dyn, dyn> interp = Gauss_lobatto(2).interpolate(basis.nodes());
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      Eigen::Map<Mat<>> qpoint_av(elems[i_elem].laplacian_av_coef(), params.n_qpoint());
      Eigen::Map<Mat<>> vert_av(&elems[i_elem].vertex_elwise_av(0), params.n_vertices());
      qpoint_av = math::hypercube_matvec(interp, vert_av);
    }
  }
}

void Solver::set_art_visc_admis() {
  stopwatch["set art visc"].stopwatch.start();
  use_art_visc = true;
  // compute the desired artificial viscosity in each element
  double char_speed = _namespace->get<double>("freestream_speed") + _namespace->get<double>("freestream_sound_speed");
  stabilizing_art_visc(_kernel_mesh(), char_speed);
  // enforce C^0 continuity
  share_vertex_data([](Element& elem, int){return elem.uncertainty;},
                    [](Element& elem, int i_vert)->double&{return elem.vertex_elwise_av(i_vert);},
                    true);
  Mat<dyn, dyn> interp = Gauss_lobatto(2).interpolate(basis.nodes());
  // interpolate from vertices to quadrature points
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Eigen::Map<Mat<>> qpoint_av(elems[i_elem].laplacian_av_coef(), params.n_qpoint());
    Eigen::Map<Mat<>> vert_av(&elems[i_elem].vertex_elwise_av(0), params.n_vertices());
    qpoint_av = math::hypercube_matvec(interp, vert_av);
  }
  stopwatch["set art visc"].stopwatch.pause();
  stopwatch["set art visc"].work_units_completed += elems.size();
}

void Solver::set_art_visc_row_size(int row_size) {
  HEXED_ASSERT(row_size >= 2, "`row_size` must be >= 2");
  HEXED_ASSERT(row_size <= basis.row_size, "`row_size` must be <= discretization row size");
  av_rs = row_size;
}

void Solver::set_fix_admissibility(bool value) {
  fix_admis = value;
}

void Solver::set_uncertainty(const Element_func& func) {
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].uncertainty = func(elems[i_elem], basis, _namespace->get<double>("flow_time"))[0];
  }
}

void Solver::set_uncertainty(double value) {
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].uncertainty = value;
  }
}

void Solver::set_uncert_surface_rep(int bc_sn) {
  const int nd = params.n_dim;
  Mat<> orth = basis.orthogonal(basis.row_size - 2).cwiseProduct(basis.node_weights());
  Mat<> weights = math::pow_outer(basis.node_weights(), nd - 1);
  auto& elems = acc_mesh->elements();
  //#pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Element& elem = elems[i_elem];
    elem.uncertainty = 0;
    Array<double> pos = elem.position(basis);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      for (int j_dim = 0; j_dim < nd; ++j_dim) {
        Mat<> ho_component = math::dimension_matvec(orth.transpose(), pos(i_dim).vector(), j_dim);
        elem.uncertainty += ho_component.dot(weights.cwiseProduct(ho_component))/(nd*nd);
      }
    }
    elem.uncertainty = std::sqrt(elem.uncertainty);
  }
}

void Solver::_update_recursive(int preti_level, double safety) {
  if (preti_level + 1 < (int)_preti_masks.size()) _update_recursive(preti_level + 1, safety);
  if (_preti_masks[preti_level]->repeat) {
    Kernel_mesh& km = _preti_masks[preti_level]->kernel_mesh;
    double dt = std::min(max_dt(safety, safety, _namespace->get<double>("time_step_limit_threshold")),
                         _namespace->get<double>("max_time_step"));
    HEXED_ASSERT(!std::isnan(dt), "time step is NaN", assert::Numerical_exception);
    bool fixed = false;
    // compute inviscid update
    for (int i = 0; i < 2; ++i) {
      Kernel_options opts {
        .sw_car = stopwatch["cartesian"],
        .sw_def = stopwatch["deformed"],
        .sw_pr = stopwatch["prolong/restrict"],
        .dt = dt,
        .i_stage = i,
        .compute_residual = false,
        .use_filter = bool(_namespace->get<int>("use_filter")),
        .mask = preti_level,
        .conv_substep = false,
      };
      apply_state_bcs();
      if (use_ldg() && !i) compute_navier_stokes(km, opts, [this](){apply_flux_bcs();}, visc, therm_cond, _namespace->get<int>("iteration")%100000 == 0 && _namespace->get<int>("iteration") != 0);
      else compute_euler(km, opts);
      // note that function call must come first to ensure it is evaluated despite short-circuiting
      fixed = fix_admissibility(_namespace->get<double>("fix_admis_max_safety"), 0) || fixed;
      stopwatch.work_units_completed += km.elems.size();
      stopwatch["cartesian"].work_units_completed += km.car_elems.size();
      stopwatch["deformed" ].work_units_completed += km.def_elems.size();
    }
    // update status for reporting
    _namespace->assign<double>("time_step", dt);
    _namespace->assign<double>("flow_time", _namespace->get<double>("flow_time") + dt);
    status.time_step = dt;
    status.flow_time += dt;
    if (preti_level + 1 < (int)_preti_masks.size()) _update_recursive(preti_level + 1, safety);
  }
}

void Solver::update() {
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  double safety = _namespace->get<double>("max_safety");
  if (_namespace->get<int>("preti")) {
    _update_recursive(0, safety);
  } else {
    double cheby_safety = _namespace->get<double>("cheby_safety");
    int inner = 0;
    for (int i_flow = 0; i_flow < _namespace->get<int>("flow_iters"); ++i_flow) {
      // compute time step
      double dt = 0;
      HEXED_ASSERT(_preti_masks.size(), "meshing mask list is empty");
      int n_preti = (_namespace->get<int>("bl_multirate") && !i_flow) ? _preti_masks.size() : 1;
      for (int i_preti = 0; i_preti < n_preti; ++i_preti) {
        int n_bl = i_preti ? _namespace->get<int>("bl_iters") : 1;
        for (int i_bl = 0; i_bl < n_bl; ++i_bl) {
          int n_cheby = i_preti ? _namespace->get<int>("n_cheby_bl") : _namespace->get<int>("n_cheby_flow");
          int max_sub_iters = i_preti ? _namespace->get<int>("max_conv_sub_iters") : 1;
          double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1, cheby_safety);
          // run chebyshev iterations
          for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby) {
            _preti_level = i_preti;
            Kernel_mesh& km = _preti_masks[_preti_level]->kernel_mesh;
            double cheby_step = math::chebyshev_step(n_cheby, i_cheby, cheby_safety);
            int sub_iters = std::ceil(max_sub_iters*cheby_step/max_cheby - 1e-6);
            double lim_thresh = n_cheby > 1 ? -1. : _namespace->get<double>("time_step_limit_threshold");
            double nominal_dt = std::min(max_dt(safety/max_cheby*sub_iters, safety, lim_thresh),
                                         _namespace->get<double>("max_time_step"));
            dt = nominal_dt*cheby_step;
            HEXED_ASSERT(!std::isnan(dt), "time step is NaN", assert::Numerical_exception);
            bool fixed = false;
            Implicit_options implicit_opts;
            if (is_implicit(_time_scheme)) {
              implicit_opts.is_implicit = true;
              implicit_opts.time_step = _namespace->get<double>("time_step");
              if (_time_scheme == crank_nicolson) implicit_opts.time_step *= .5;
              if (_time_scheme == dirk2) implicit_opts.time_step *= dirk2_gamma;
            }
            for (int i_sub = 0; i_sub < sub_iters; ++i_sub) {
              // compute inviscid update
              for (int i = 0; i < 2; ++i) {
                Kernel_options opts {
                  .sw_car = stopwatch["cartesian"],
                  .sw_def = stopwatch["deformed"],
                  .sw_pr = stopwatch["prolong/restrict"],
                  .dt = dt/sub_iters,
                  .i_stage = i,
                  .compute_residual = false,
                  .use_filter = bool(_namespace->get<int>("use_filter")),
                  .mask = i_preti,
                  .conv_substep = (sub_iters > 1) && use_ldg(),
                  .implicit_opts = implicit_opts,
                };
                apply_state_bcs();
                if (use_ldg() && !i && !i_sub) {
                  compute_navier_stokes(km, opts, [this](){apply_flux_bcs();}, visc, therm_cond, true);
                } else {
                  compute_euler(km, opts);
                }
                ++inner;
                // note that function call must come first to ensure it is evaluated despite short-circuiting
                bool f = fix_admissibility(_namespace->get<double>("fix_admis_max_safety"), inner);
                fixed = f || fixed;
                if (f) break;
              }
              stopwatch.work_units_completed += km.elems.size();
              stopwatch["cartesian"].work_units_completed += km.car_elems.size();
              stopwatch["deformed" ].work_units_completed += km.def_elems.size();
            }
            // update status for reporting
            if (!is_implicit(_time_scheme)) {
              _namespace->assign<double>("time_step", dt);
              _namespace->assign<double>("flow_time", _namespace->get<double>("flow_time") + dt);
              status.time_step = dt;
              status.flow_time += dt;
            }
          }
        }
      }
    }
  }
  ++status.iteration;
  stopwatch.stopwatch.pause();
}

void Solver::smooth_init_cond(Int n_iter) {
  printers::info("Smoothing initial state...");
  auto km = _kernel_mesh();
  Kernel_options opts {
    stopwatch["fix admis."]["cartesian"],
    stopwatch["fix admis."]["deformed"],
    stopwatch["prolong/restrict"],
    1.,
    0,
  };
  max_dt_fix_therm_admis(km, opts, 1., _namespace->get<double>("fix_admis_max_safety"), true);
  for (Int iter = 0; iter < n_iter; ++iter) {
    apply_state_bcs();
    compute_fix_therm_admis(km, opts, [this](){apply_flux_bcs();});
  }
  printers::info(" done.\n");
}

void Solver::compute_residual() {
  #if 0
  apply_state_bcs();
  auto compute_discon = [this](bool is_flux) {
    Mat<> face_weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
    int nd = params.n_dim;
    Array<double> face_min = Array<double>::make_uniform({params.n_var}, huge);
    Array<double> face_max = Array<double>::make_uniform({params.n_var}, -huge);
    for (bool is_def : {0, 1}) {
      #pragma omp parallel for reduction(min:face_min) reduction(max:face_max)
      for (auto& con : acc_mesh->neighbor_connections(is_def)) {
        for (int i_side = 0; i_side < 2; ++i_side) {
          Array<double> state = con.face(i_side).flow_state()(is_flux);
          for (int i_var = 0; i_var < params.n_var; ++i_var) {
            for (int i_qpoint = 0; i_qpoint < params.n_face_qpoint(); ++i_qpoint) {
              face_min[i_var] = std::min(face_min[i_var], state(i_var)[i_qpoint]);
              face_max[i_var] = std::max(face_max[i_var], state(i_var)[i_qpoint]);
            }
          }
        }
      }
    }
    for (bool is_def : {0, 1}) {
      #pragma omp parallel for
      for (auto& con : acc_mesh->neighbor_connections(is_def)) {
        bool farfield = false;
        for (int i_side = 0; i_side < 2; ++i_side) {
          if (con.face(i_side).boundary_connection()) {
            farfield = farfield || con.face(i_side).boundary_connection()->boundary_condition() < 2*nd;
          }
        }
        if (farfield) {
          for (int i_side = 0; i_side < 2; ++i_side) con.face(i_side).discontinuity()(is_flux) = 0;
        } else {
          Array<double> diff = con.face(1).flow_state()(is_flux).copy();
          auto dir = con.get_direction();
          if (is_flux) diff *= math::sign(dir.flip_normal(0) == dir.flip_normal(1));
          auto perm = face_permutation(nd, params.row_size, dir, diff.data(), turb);
          perm->match_faces();
          diff -= con.face(0).flow_state()(is_flux);
          for (int i_var = 0; i_var < params.n_var; ++i_var) if (i_var != params.n_dim || !is_flux) {
            double norm = std::sqrt(diff(i_var).vector().dot(face_weights.cwiseProduct(diff(i_var).vector())));
            norm /= face_max[i_var] - face_min[i_var] + 1e-15*(std::abs(face_max[i_var]) + std::abs(face_min[i_var]));
            for (int i_side = 0; i_side < 2; ++i_side) con.face(i_side).discontinuity()(is_flux)[i_var] = norm;
          }
        }
      }
    }
  };
  compute_discon(false);
  Kernel_options opts {
    .sw_car = stopwatch["cartesian"],
    .sw_def = stopwatch["deformed"],
    .sw_pr = stopwatch["prolong/restrict"],
    .dt = 1.,
    .i_stage = 0,
    .compute_residual = true,
    .use_filter = bool(_namespace->get<int>("use_filter")),
  };
  if (use_ldg()) {
    auto bc_fun = [this, compute_discon]() {
      apply_flux_bcs();
      compute_discon(true);
    };
    compute_navier_stokes(_kernel_mesh(), opts, bc_fun, visc, therm_cond, false);
  } else {
    compute_euler(_kernel_mesh(), opts);
  }
  #pragma omp parallel for
  for (auto& vec : acc_mesh->face_refinements()) {
    for (int i_ref = vec.size() - 1; i_ref >= 0; --i_ref) {
      auto& ref = vec[i_ref];
      ref.coarse().discontinuity() = .5*(ref.fine()[0]->discontinuity() + ref.fine()[1]->discontinuity());
    }
  }
  auto& elems = acc_mesh->elements();
  Mat<> weights_1d = basis.node_weights();
  Mat<> weights = math::pow_outer(weights_1d, params.n_dim);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* res = elems[i_elem].residual_cache();
    elems[i_elem].residual = 0;
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      double mean_sq = 0;
      for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
        double r = res[i_var*params.n_qpoint() + i_qpoint];
        mean_sq += r*r*weights(i_qpoint);
      }
      elems[i_elem].residual += std::sqrt(mean_sq);
    }
  }
  double cumulative_max = 0;
  int min_level = 5;
  for (int level = 0; level < (int)_preti_masks.size(); ++level) {
    double max_res = 0;
    auto km = _preti_masks[level]->kernel_mesh;
    #pragma omp parallel for reduction(max:max_res)
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      if (elems[i_elem].mask() == level) max_res = std::max(max_res, elems[i_elem].residual);
    }
    _preti_masks[level]->max_residual = max_res;
    cumulative_max = std::max(cumulative_max, max_res);
    _preti_masks[level]->repeat = !level;
  }
  for (int level = min_level; level < (int)_preti_masks.size(); ++level) {
    int desired = 10*_preti_masks[level]->max_residual/cumulative_max;
    _preti_masks[level]->desired_iters = desired;
    for (int add_at = level; add_at >= min_level && _effective_preti_iters(level) < desired; --add_at) {
      _preti_masks[add_at]->repeat = true;
    }
  }
  #else
  Kernel_options opts {
    .sw_car = stopwatch["cartesian"],
    .sw_def = stopwatch["deformed"],
    .sw_pr = stopwatch["prolong/restrict"],
    .dt = 1.,
    .i_stage = 0,
    .compute_residual = true,
    .use_filter = false,
  };
  if (use_ldg()) {
    compute_navier_stokes(_kernel_mesh(), opts, [this](){apply_flux_bcs();}, visc, therm_cond, false);
  } else {
    compute_euler(_kernel_mesh(), opts);
  }
  #endif
  Int n_diff = 0;
  Int n_source = 0;
  #pragma omp parallel for reduction(+:n_diff,n_source)
  for (Int i_elem = 0; i_elem < _preti_masks[0]->kernel_mesh.elems.size(); ++i_elem) {
    n_diff += _preti_masks[0]->kernel_mesh.elems[i_elem].diffusion_limited;
    n_source += _preti_masks[0]->kernel_mesh.elems[i_elem].source_limited;
  }
  _namespace->assign<int>("n_diffusion_limited", n_diff);
  _namespace->assign<int>("n_source_limited", n_source);
}

Int Solver::_effective_preti_iters(int level) {
  Int iters = 1;
  for (int l = 1; l <= level; ++l) {
    if (_preti_masks[l]->repeat) iters *= 2;
  }
  return iters;
}

void Solver::print_preti_iters() {
  if (!_namespace->get<int>("preti")) return;
  printers::info("PRETI sub-iterations (" + to_string((Int)_preti_masks.size()) + " levels total):\n");
  std::string line0 = "repeat?:                   ";
  std::string line1 = "total effective iterations:";
  std::string line2 = "desired iterations:        ";
  std::string line3 = "level max residual:        ";
  for (Int i_preti = 0; i_preti < (Int)_preti_masks.size(); ++i_preti) {
    line0 += format_str(" %8li", (Int)_preti_masks[i_preti]->repeat);
    line1 += format_str(" %8li", _effective_preti_iters(i_preti));
    line2 += format_str(" %8li", _preti_masks[i_preti]->desired_iters);
    line3 += format_str(" %8.1e", _preti_masks[i_preti]->max_residual);
  }
  printers::info(line0 + "\n" + line1 + "\n" + line2 + "\n" + line3 + "\n");
}


void Solver::compute_spectral_uncertainty() {
  std::vector<int> vars;
  for (int i_var = 0; i_var < params.n_dim + 2; ++i_var) vars.push_back(i_var);
  if (use_art_visc) vars.push_back(params.n_var + 3);
  int nv = vars.size();
  Array<double> state_min = Array<double>::make_uniform({nv}, huge);
  Array<double> state_max = Array<double>::make_uniform({nv}, -huge);
  auto& elems = acc_mesh->elements();
  #pragma omp parallel for reduction(min:state_min) reduction(max:state_max)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Array<double> state = elems[i_elem].numeric_state();
    for (int i_var = 0; i_var < nv; ++i_var) {
      state_min[i_var] = std::min(state_min[i_var], state(vars[i_var]).extreme(0));
      state_max[i_var] = std::max(state_min[i_var], state(vars[i_var]).extreme(1));
    }
  }
  Mat<dyn, dyn> orth = basis.orthogonal(params.row_size - 1).transpose()*basis.node_weights().asDiagonal();
  Mat<> weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
  #pragma omp parallel
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    Array<double> state = elem.numeric_state();
    elem.spectral_uncert() = 0;
    elem.flux_uncert = 0;
    for (int i_var = 0; i_var < nv; ++i_var) {
      for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
        Mat<> proj = math::dimension_matvec(orth, state(vars[i_var]).vector(), i_dim);
        double normalize = state_max[i_var] - state_min[i_var];
        if (i_var >= params.n_var) normalize *= 10;
        double& elem_uncert = elem.spectral_uncert()[i_dim];
        elem_uncert = std::max(elem_uncert, std::sqrt(proj.dot(proj.cwiseProduct(weights)))/normalize);
      }
    }
  }
  double total_sq_flux = 1;
  double total_area = 1;
  if (visc.is_viscous) {
    total_sq_flux = 0;
    total_area = 0;
    auto bc_fun = [this, &total_sq_flux, &total_area]() {
      apply_flux_bcs();
      Array<double> face_weights(math::pow_outer(basis.node_weights(), params.n_dim - 1));
      #pragma omp parallel for reduction(+:total_sq_flux,total_area)
      for (auto& con : acc_mesh->neighbor_connections(true)) {
        auto dir = con.get_direction();
        if (!con.has_elements() || dir.i_dim[0] != dir.i_dim[1]) continue;
        if (!con.face(0).element()->is_extruded() || !con.face(1).element()->is_extruded()) continue;
        if (dir.i_dim[0] != con.face(0).element()->wall_dimension()) continue;
        Array<double> flux_diff = con.face(0).flow_state()(1) - con.face(1).flow_state()(1);
        Array<double> flux_avg  = con.face(0).flow_state()(1) + con.face(1).flow_state()(1);
        Array<double> nrml = con.face(0).normal();
        Array<double> area = Array<double>::make_uniform({params.n_face_qpoint()}, 0.);
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          area += nrml(i_dim)*nrml(i_dim);
        }
        area.sqrt(true);
        area *= con.face(0).nominal_area();
        double uncert = 0;
        double total = 0;
        for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
          uncert += (flux_diff(i_dim)*flux_diff(i_dim)*face_weights/(area*area)).sum();
          total += (flux_avg(i_dim)*flux_avg(i_dim)*face_weights/area).sum();
        }
        for (int i_side = 0; i_side < 2; ++i_side) {
          con.face(i_side).element()->flux_uncert += std::sqrt(uncert);
        }
        if (con.face(0).element()->has_wall() || con.face(1).element()->has_wall()) {
          total_sq_flux += total;
          total_area += (area*face_weights).sum();
        }
      }
    };
    Kernel_options opts {
      .sw_car = stopwatch["cartesian"],
      .sw_def = stopwatch["deformed"],
      .sw_pr = stopwatch["prolong/restrict"],
      .dt = 1.,
      .i_stage = 0,
      .compute_residual = true,
      .use_filter = bool(_namespace->get<int>("use_filter")),
    };
    compute_navier_stokes(_kernel_mesh(), opts, bc_fun, visc, therm_cond, false);
  }
  _namespace->assign("rms_flux", std::sqrt(total_sq_flux/total_area));
}

void Solver::update_bound_conds() {
  double relative = _namespace->get<double>("max_roughness_relative")*_namespace->get<double>("geom_length");
  double max_rough = std::min(_namespace->get<double>("max_roughness_absolute"), relative);
  if (_namespace->get<int>("local_roughness")) {
    _namespace->assign("hexed_max_roughness", max_rough);
  } else {
    std::string expr = "temperature = energy*(heat_rat - 1)/specific_gas_air; $transport_expr;"
                       "inv_roughness = sqrt(sqrt(visc_stress0^2 + visc_stress1^2 + visc_stress2^2)/density)*density"
                       "/(dynamic_viscosity*max_roughness_plus);";
    bounds_surface(expr, 2*params.n_dim, 20);
    if (_namespace->exists("max_surface_inv_roughness")) {
      double rough = std::min(1./_namespace->get<double>("max_surface_inv_roughness"), max_rough);
      _namespace->assign("hexed_surface_roughness", rough);
    }
  }
  auto inter {_interpreter()};
  auto bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    int bc_sn = con.boundary_condition();
    acc_mesh->boundary_condition(bc_sn).set_prescribed(inter, con);
  }
}

void Solver::compute_lts_constraints() {
  auto& elems = acc_mesh->deformed().elements();
  int nd = params.n_dim;
  int nq = params.n_qpoint();
  int n_cheby = _namespace->get<int>("n_cheby_bl");
  double cheby_safety = _namespace->get<double>("cheby_safety");
  double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1, cheby_safety);
  // write local time steps for convection and diffusion to the mass and energy of the reference state.
  // Reference state is used for storage because `Element::time_step_scale` only has space for one scalar
  for (int i_term = 0; i_term < 2; ++i_term) {
    double safeties [] {1., huge};
    max_dt(safeties[i_term]/max_cheby, safeties[!i_term], -1.);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      Eigen::Map<Mat<>>(elems[i_elem].residual_cache() + (params.n_dim + i_term)*nq, nq) = Eigen::Map<Mat<>>(elems[i_elem].time_step_scale(), nq);
    }
  }
  double min_ratio = huge;
  #pragma omp parallel for reduction(min:min_ratio)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* cache = elems[i_elem].residual_cache();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      min_ratio = std::min(min_ratio, cache[(nd + 1)*nq + i_qpoint]/cache[nd*nq + i_qpoint]);
    }
  }
  _namespace->assign("min_lts_dc_ratio", min_ratio);
}

Iteration_status Solver::iteration_status() {
  Iteration_status stat = status;
  return stat;
}

bool Solver::is_admissible() {
  auto& sw = stopwatch["fix admis."]["check admis."];
  sw.stopwatch.start();
  auto& elems = _preti_masks[_preti_level]->kernel_mesh.elems;
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  bool admiss = true;
  bool finite = true;
  std::string message;
  auto check_admis = [&](double* data, int n_qpoint, int n_var) {
    bool adm = true;
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      adm = adm && (data[nd*n_qpoint + i_qpoint] > 0.)
                && (data[(nd + 1)*n_qpoint + i_qpoint] > 0.);
      for (int i_var = 0; i_var < n_var; ++i_var) {
        if (!std::isfinite(data[i_var*n_qpoint + i_qpoint])) {
          finite = false;
          #pragma omp critical
          message = format_str("variable %i = %e has non-finite value.", i_var, data[i_var*n_qpoint + i_qpoint]);
        }
      }
    }
    return adm;
  };
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].record = 0;
  }
  #pragma omp parallel for reduction(&&:admiss,finite)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    bool elem_admis = true;
    elem_admis = elem_admis && check_admis(elem.state(), nq, params.n_var);
    for (int i_face = 0; i_face < params.n_dim*2; ++i_face) {
      elem_admis = elem_admis && check_admis(elem.face(i_face, false), nq/rs, nd + 2);
    }
    if (!elem_admis) elem.record = 1;
    admiss = admiss && elem_admis;
  }
  auto& face_refs = _preti_masks[_preti_level]->kernel_mesh.face_refinements;
  bool refined_admiss = 1;
  #pragma omp parallel for reduction (&&:refined_admiss,finite)
  for (auto& vec : face_refs) {
    for (auto& ref : vec) {
      for (int i_fine = 0; i_fine < 2; ++i_fine) {
        refined_admiss = refined_admiss && check_admis(ref.fine[i_fine][0], nq/rs, nd + 2);
      }
    }
  }
  HEXED_ASSERT(finite, message, assert::Numerical_exception)
  sw.work_units_completed += acc_mesh->elements().size();
  sw.stopwatch.pause();
  return admiss && refined_admiss;
}

bool Solver::fix_admissibility(double stability_ratio, int sub_iter) {
  if (!fix_admis) return false;
  auto& sw_fix = stopwatch["fix admis."];
  sw_fix.stopwatch.start();
  std::string wd = _namespace->get<std::string>("working_dir");
  std::string vis_expr = _namespace->get<std::string>("vis_field_vars");
  std::vector<double> freestream(params.n_var);
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    freestream[i_var] = _namespace->get<double>(str_cat("freestream", i_var));
  }
  std::vector<double> sanity_interval(params.n_var, huge);
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    double stag_enthalpy = freestream[params.n_dim + 1] + _namespace->get<double>("freestream_pressure");
    sanity_interval[i_dim] = 1e3*std::sqrt(freestream[params.n_dim]*stag_enthalpy);
  }
  for (int i_var : {params.n_dim, params.n_dim + 1}) sanity_interval[i_var] = 1e3*freestream[i_var];
  if (turb == k_omega) {
    // TKE should be less than total energy
    sanity_interval[params.n_dim + 2] = freestream[params.n_dim + 1];
    // the specific turbulent dissipation varying by a factor of more than exp(1e3) would be implausible
    sanity_interval[params.n_dim + 3] = freestream[params.n_dim]*1e3;
  }
  int nq = params.n_qpoint();
  auto& elems = _preti_masks[_preti_level]->kernel_mesh.elems;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].state();
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        double& s = state[i_var*nq + i_qpoint];
        if (!std::isfinite(s)) {
          s = freestream[i_var];
        } else if (s > freestream[i_var] + sanity_interval[i_var]) {
          s = freestream[i_var] + sanity_interval[i_var];
        } else if (s < freestream[i_var] - sanity_interval[i_var]) {
          s = freestream[i_var] - sanity_interval[i_var];
        }
      }
    }
  }
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  int iter = 0;
  int n_iters = std::numeric_limits<int>::max();
  for (; iter < n_iters;) {
    HEXED_ASSERT(iter < 100'000, format_str("failed to fix thermodynamic admissability in %i iterations", iter),
                 assert::Numerical_exception)
    if (is_admissible()) {
      if (iter) {
        n_iters = std::min(n_iters, 2*iter);
      } else {
        break;
      }
    } else {
      n_iters = std::numeric_limits<int>::max();
    }
    for (int i_vis = 0; i_vis < 2; ++i_vis) {
      if (iter == (i_vis + 1)*1000) {
        double ft = _namespace->get<double>("flow_time");
        _namespace->assign<double>("flow_time", i_vis);
        visualize_field("default", str_cat(wd, "severe_inadmis", status.iteration, "_", sub_iter, "_", i_vis),
                        vis_expr);
        _namespace->assign("flow_time", ft);
      }
    }
    if (iter == 0) {
      printers::warn("Warning: ", true);
      printers::warn(str_cat("Nonphysical flow state detected (solver iteration ", status.iteration,
                             " sub-iteration ", sub_iter, "). Attempting to fix...\n"));
    }
    printers::warn(format_str("    iteration %i\n", iter));
    if (status.iteration >= last_fix_vis_iter + 1000 && iter == 0) {
      last_fix_vis_iter = status.iteration;
      visualize_field("default", str_cat(wd, "inadmis", status.iteration, "_", sub_iter), vis_expr);
    }
    for (int inner = 0; inner < 100; ++inner, ++iter) {
      double dt = stability_ratio;
      Kernel_options opts {
        stopwatch["fix admis."]["cartesian"],
        stopwatch["fix admis."]["deformed"],
        stopwatch["prolong/restrict"],
        0.,
        0,
        false,
        false,
      };
      max_dt_fix_therm_admis(_kernel_mesh(), opts, dt, dt, true);
      #if 1
      auto bc_cons {acc_mesh->boundary_connections()};
      #pragma omp parallel for
      for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
        bc_cons[i_con].ghost().flow_state()(0) = bc_cons[i_con].inside().flow_state()(0);
      }
      opts.dt = 1.;
      compute_fix_therm_admis(_kernel_mesh(), opts, [this](){apply_fta_flux_bcs();});
      #else
      apply_state_bcs();
      opts.dt = 1.;
      compute_fix_therm_admis(_kernel_mesh(), opts, [this](){apply_flux_bcs();});
      #endif
    }
  }
  if (iter) printers::warn("done\n");
  status.fix_admis_iters += iter;
  _namespace->assign("fix_iters", _namespace->get<int>("fix_iters") + iter);
  sw_fix.work_units_completed += acc_mesh->elements().size()*iter;
  sw_fix.stopwatch.pause();
  return iter;
}

void Solver::reset_counters() {
  status.fix_admis_iters = 0;
  _namespace->assign("fix_iters", 0);
}

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint,
                                   const Qpoint_func& func) {
  return func(acc_mesh->element(ref_level, is_deformed, serial_n), basis, i_qpoint,
                                _namespace->get<double>("flow_time"));
}

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, const Element_func& func) {
  return func(acc_mesh->element(ref_level, is_deformed, serial_n), basis, _namespace->get<double>("flow_time"));
}

std::vector<double> Solver::integral_field(const Qpoint_func& integrand) {
  Stopwatch_tree::Starter sw_starter(stopwatch["integrals"]["field"]);
  // compute `n_dim`-dimensional quadrature weights from 1D weights
  Eigen::VectorXd weights = math::pow_outer(basis.node_weights(), params.n_dim);
  // now compute the integral with the above quadrature weights
  Mat<dyn, dyn> integral = Mat<dyn, dyn>::Zero(integrand.n_var(params.n_dim), 1);
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for reduction(+:integral)
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    Element& element {elements[i_elem]};
    double volume = element.nominal_volume();
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
      auto qpoint_integrand {integrand(element, basis, i_qpoint, _namespace->get<double>("flow_time"))};
      for (unsigned i_var = 0; i_var < qpoint_integrand.size(); ++i_var) {
        integral(i_var) += weights[i_qpoint]*volume*qpoint_integrand[i_var]*element.jacobian_determinant(i_qpoint);
      }
    }
  }
  stopwatch["integrals"]["field"].work_units_completed += 1;
  stopwatch["integrals"].work_units_completed += 1;
  return {integral.data(), integral.data() + integral.size()};
}

//! \cond
template <typename T>
class Vis_evaluator {
  public:
  Vis_evaluator(Interpreter&& inter, std::function<void(Namespace&, T&)> assign, std::string expr, T& t,
                int n_dim_topo)
  : _inter{inter}, _assign{assign}, _expr{expr}, _n_dim_topo{n_dim_topo}
  {
    Storage_params params = t.storage_params();
    _n_dim = params.n_dim;
    auto sub = _inter.make_sub();
    _assign(*sub.variables, t);
    sub.subspace();
    sub.exec(_expr);
    _var_names = sub.variables->names();
    std::erase_if(_var_names, [&sub](std::string name) {
      return !sub.variables->lookup<double>(name) && !sub.variables->lookup<Array<double>>(name);
    });
    _n_var = _var_names.size();
    _shape = hypercubes(_n_dim + _n_var, _n_dim_topo, params.row_size);
  }

  Array<double> evaluate(T& t) {
    Array<double> qpoints(_shape);
    auto sub = _inter.make_sub();
    _assign(*sub.variables, t);
    sub.exec(_expr);
    for (int i_dim = 0; i_dim < _n_dim; ++i_dim) {
      qpoints(i_dim) = sub.variables->template get<Array<double>>("pos" + std::to_string(i_dim));
    }
    for (int i_var = 0; i_var < _n_var; ++i_var) {
      sub.variables->assign_array(qpoints(_n_dim + i_var), _var_names[i_var]);
    }
    return qpoints;
  }
  std::vector<std::string> var_names() {return _var_names;}

  void visualize(std::string format, std::string name, int n_sample, bool wireframe, next::Sequence<T&> seq,
                 double time, const Basis& basis, std::function<bool(T&)> mask) {
    auto visualizer = Visualizer::create(format, _n_dim, wireframe ? 1 : _n_dim_topo, name, _var_names,
                                         time, Visualizer::block);
    #pragma omp parallel for
    for (Int i = 0; i < seq.size(); ++i) if (mask(seq[i])) {
      Array<double> qpoints {evaluate(seq[i])};
      Vis_data vis_dat(qpoints, basis);
      if (wireframe) {
        Array<double> edges {vis_dat.edges(n_sample)};
        for (int i_dim = 0; i_dim < _n_dim_topo; ++i_dim) {
          for (int i_edge = 0; i_edge < edges(i_dim).shape()[0]; ++i_edge) {
            #pragma omp critical
            visualizer->write_block(edges(i_dim)(i_edge)(0, _n_dim),
                                    edges(i_dim)(i_edge)(_n_dim, _n_dim + _n_var));
          }
        }
      } else {
        Array<double> interior {vis_dat.interior(n_sample)};
        #pragma omp critical
        visualizer->write_block(interior(0, _n_dim), interior(_n_dim, _n_dim + _n_var));
      }
    }
  }

  private:
  Interpreter _inter;
  std::function<void(Namespace&, T&)> _assign;
  std::string _expr;
  Int _n_var;
  Int _n_dim;
  Int _n_dim_topo;
  std::vector<Int> _shape;
  std::vector<std::string> _var_names;
};
//! \endcond

void Solver::bounds_surface(std::string expr, int bc_sn, int n_sample = 20) {
  // setup
  const int nd = params.n_dim;
  auto bc_cons {acc_mesh->boundary_connections()};
  Boundary_connection* reference_con = nullptr;
  for (auto& con : bc_cons) {
    if (con.boundary_condition() == bc_sn) reference_con = &con;
  }
  if (!reference_con) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, *reference_con, params.n_dim - 1
  );
  std::vector<std::string> var_names = evaluator.var_names();
  Int n_var = var_names.size();
  Mat<> weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
  // write the state to the faces so that the BCs can access it
  compute_write_face(_kernel_mesh());
  // compute the integral
  Array<double> bounds({2, n_var});
  for (int i_var = 0; i_var < n_var; ++i_var) {
    bounds(0)[i_var] = huge;
    bounds(1)[i_var] = -huge;
  }
  //#pragma omp parallel for reduction(+:integral)
  for (Int i_con = 0; i_con < (Int)bc_cons.size(); ++i_con) {
    auto& con {bc_cons[i_con]};
    if (con.boundary_condition() != bc_sn) continue;
    Array<double> qpoints{evaluator.evaluate(con)};
    Vis_data vis_dat(qpoints(nd, end), basis);
    Array<double> interior = vis_dat.interior(n_sample);
    for (int i_var = 0; i_var < n_var; ++i_var) {
      for (int i_point = 0; i_point < interior(i_var).size(); ++i_point) {
        bounds(0)[i_var] = std::min(bounds(0)[i_var], interior(i_var)[i_point]);
        bounds(1)[i_var] = std::max(bounds(1)[i_var], interior(i_var)[i_point]);
      }
    }
  }
  for (int i_var = 0; i_var < n_var; ++i_var) {
    _namespace->assign("min_surface_" + var_names[i_var], bounds(0)[i_var]);
    _namespace->assign("max_surface_" + var_names[i_var], bounds(1)[i_var]);
  }
}

void Solver::integrate_field(std::string expr) {
  Stopwatch_tree::Starter sw_starter(stopwatch["integrals"]["field"]);
  // setup
  const int nd = params.n_dim;
  auto& elems {acc_mesh->elements()};
  if (!elems.size()) return;
  Vis_evaluator<Element> evaluator(
    _interpreter(),
    [&](Namespace& space, Element& elem) {vis_variables::field(space, elem, basis);},
    expr, elems[0], params.n_dim
  );
  std::vector<std::string> var_names = evaluator.var_names();
  Mat<> weights = math::pow_outer(basis.node_weights(), params.n_dim);
  // compute the integral
  Mat<dyn, dyn> integral = Mat<dyn, dyn>::Zero(var_names.size(), 1);
  #pragma omp parallel for reduction(+:integral)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    double volume = math::pow(elem.nominal_size(), params.n_dim);
    Array<double> qpoints{evaluator.evaluate(elem)};
    for (int i_var = 0; i_var < (int)var_names.size(); ++i_var) {
      integral(i_var) += volume*weights.dot(qpoints(nd + i_var).vector());
    }
  }
  for (int i_var = 0; i_var < (int)var_names.size(); ++i_var) {
    _namespace->assign("integral_field_" + var_names[i_var], integral(i_var));
  }
  ++stopwatch["integrals"]["field"].work_units_completed;
  ++stopwatch["integrals"].work_units_completed;
};

void Solver::integrate_surface(std::string expr, int bc_sn) {
  Stopwatch_tree::Starter sw_starter(stopwatch["integrals"]["surface"]);
  // setup
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int nfq = nq/basis.row_size;
  auto bc_cons {acc_mesh->boundary_connections()};
  Boundary_connection* reference_con = nullptr;
  for (auto& con : bc_cons) {
    if (con.boundary_condition() == bc_sn) reference_con = &con;
  }
  if (!reference_con) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, *reference_con, params.n_dim - 1
  );
  std::vector<std::string> var_names = evaluator.var_names();
  Mat<> weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
  // write the state to the faces so that the BCs can access it
  compute_write_face(_kernel_mesh());
  // compute the integral
  Mat<dyn, dyn> integral = Mat<dyn, dyn>::Zero(var_names.size(), 1);
  #pragma omp parallel for reduction(+:integral)
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con {bc_cons[i_con]};
    if (con.boundary_condition() != bc_sn) continue;
    double area = con.inside().nominal_area();
    Array<double> qpoints{evaluator.evaluate(con)};
    Array<double> nrml = con.normal().copy();
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      double nrml_mag = 0;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        nrml_mag += math::pow(nrml(i_dim)[i_qpoint], 2);
      }
      nrml_mag = std::sqrt(nrml_mag);
      for (int i_var = 0; i_var < (int)var_names.size(); ++i_var) {
        integral(i_var) += nrml_mag*weights(i_qpoint)*area*qpoints(nd + i_var)[i_qpoint];
      }
    }
  }
  for (int i_var = 0; i_var < (int)var_names.size(); ++i_var) {
    _namespace->assign("integral_surface_" + var_names[i_var], integral(i_var));
  }
  ++stopwatch["integrals"]["surface"].work_units_completed;
  ++stopwatch["integrals"].work_units_completed;
};

void Solver::visualize_field(std::string format, std::string name, std::string expr, int n_sample, bool wireframe) {
  std::string sw_name = "field";
  if (wireframe) sw_name = sw_name + " wireframe";
  Stopwatch_tree::Starter sw_starter(stopwatch["visualization"][sw_name]);
  HEXED_ASSERT(params.n_dim > wireframe, "can only visualize field wireframes in > 1D");
  auto& elems = acc_mesh->elements();
  if (!elems.size()) return;
  Vis_evaluator<Element> evaluator(
    _interpreter(),
    [&](Namespace& space, Element& elem) {vis_variables::field(space, elem, basis);},
    expr, elems[0], params.n_dim
  );
  next::Sequence<Element&> elem_seq = {
    [&elems](Int index)->Element& {return elems[index];},
    [&elems]()->Int {return elems.size();},
  };
  evaluator.visualize(format, name, n_sample, wireframe, elem_seq,
                      _namespace->get<double>("flow_time"), basis, [](Element&){return true;});
  ++stopwatch["visualization"].work_units_completed;
  stopwatch["visualization"][sw_name].work_units_completed += elems.size();
}

void Solver::visualize_surface(std::string format, std::string name, int bc_sn, std::string expr,
                               int n_sample, bool wireframe) {
  std::string sw_name = "surface";
  if (wireframe) sw_name = sw_name + " wireframe";
  Stopwatch_tree::Starter sw_starter(stopwatch["visualization"][sw_name]);
  HEXED_ASSERT(params.n_dim > 1, "cannot visualize surfaces in 1D");
  HEXED_ASSERT(params.n_dim > 1 + wireframe, "can only visualize surface wireframes in 3D");
  auto bc_cons {acc_mesh->boundary_connections()};
  Boundary_connection* reference_con = nullptr;
  for (auto& con : bc_cons) {
    if (con.boundary_condition() == bc_sn) reference_con = &con;
  }
  if (!reference_con) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, *reference_con, params.n_dim - 1
  );
  evaluator.visualize(format, name, n_sample, wireframe, bc_cons,
                      _namespace->get<double>("flow_time"), basis,
                      [bc_sn](Boundary_connection& con){return con.boundary_condition() == bc_sn;});
  ++stopwatch["visualization"].work_units_completed;
  stopwatch["visualization"][sw_name].work_units_completed += bc_cons.size();
}

void Solver::visualize_contour(std::string format, std::string name, std::string contour_expr,
                               std::string vis_expr, double const_tol, int n_sample) {
  Stopwatch_tree::Starter sw_starter(stopwatch["visualization"]["contour"]);
  HEXED_ASSERT(params.n_dim > 1, "cannot visualize surfaces in 1D");
  auto& elems = acc_mesh->elements();
  if (!elems.size()) return;
  vis_expr = vis_expr + ";hexed_contour = " + contour_expr + ";";
  Vis_evaluator<Element> evaluator(
    _interpreter(),
    [&](Namespace& space, Element& elem) {vis_variables::field(space, elem, basis);},
    vis_expr, elems[0], params.n_dim
  );
  auto var_names = evaluator.var_names();
  int i_contour = std::find(var_names.begin(), var_names.end(), "hexed_contour") - var_names.begin();
  auto visualizer = Visualizer::create(format, params.n_dim, params.n_dim - 1, name, var_names,
                                       _namespace->get<double>("flow_time"), Visualizer::block);
  Int n_write = 0;
  #pragma omp parallel for reduction(+:n_write)
  for (Int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Array<double> qpoints {evaluator.evaluate(elems[i_elem])};
    Vis_data data(qpoints, basis);
    auto contour = data.compute_contour(params.n_dim + i_contour, 0., n_sample/2, 4, const_tol);
    if (contour.elem_vert_inds.size()) {
      ++n_write;
      Array<double> values = data.sample(contour.vert_ref_coords);
      #pragma omp critical
      visualizer->write_unstruct(contour.elem_vert_inds, values(0, params.n_dim), values(params.n_dim, end));
    }
  }
  if (!n_write) {
    printers::warn("Warning: ", true);
    printers::warn("Contour is empty. You won't be able to open it in Paraview. ");
  }
  ++stopwatch["visualization"].work_units_completed;
  stopwatch["visualization"]["contour"].work_units_completed += n_write;
}

void Solver::vis_lts_constraints(std::string format, std::string name, int n_sample) {
  auto& elems = acc_mesh->elements();
  int nf = params.n_dof();
  int nq = params.n_qpoint();
  // write local time steps for convection and diffusion to the mass and energy of the reference state.
  // Reference state is used for storage because `Element::time_step_scale` only has space for one scalar
  for (int i_term = 0; i_term < 2; ++i_term) {
    double safeties [] {1., huge};
    max_dt(safeties[i_term], safeties[!i_term], -1.);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      Eigen::Map<Mat<>>(elems[i_elem].residual_cache() + (params.n_dim + i_term)*nq, nq) = Eigen::Map<Mat<>>(elems[i_elem].time_step_scale(), nq);
    }
  }
  // swap current state and reference state
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Eigen::Map<Mat<>> state(elems[i_elem].state(), nf);
    Eigen::Map<Mat<>> res_cache(elems[i_elem].residual_cache(), nf);
    Mat<> temp = res_cache;
    res_cache = state;
    state = temp;
  }
  // visualize. Note that visualizing straight from the reference state would require implementing another `Qpoint_func` which would be ugly
  std::string expr {"lts_convective = density; lts_diffusive = energy; lts_ratio = lts_diffusive/lts_convective;"};
  visualize_field(format, name, expr, n_sample);
  // restore the current state from the reference state
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Eigen::Map<Mat<>>(elems[i_elem].state(), nf) = Eigen::Map<Mat<>>(elems[i_elem].residual_cache(), nf);
  }
}

Array<double> Solver::skews() {
  auto& elems = acc_mesh->elements();
  Array<double> s({elems.size()});
  Equiangle_skewness equi;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    s[i_elem] = equi(elems[i_elem], basis, _namespace->get<double>("flow_time"))[0];
  }
  return s;
}

}
