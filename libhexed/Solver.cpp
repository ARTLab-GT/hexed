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

void Solver::share_vertex_data(std::function<double&(Element&, int i_vertex)> access_fun,
                               Solver::Reduction reduction) {
  share_vertex_data(access_fun, access_fun, reduction);
}

void Solver::share_vertex_data(std::function<double(Element&, int i_vertex)> get,
                               std::function<double&(Element&, int i_vertex)> set,
                               Solver::Reduction reduction) {
  int nv = params.n_vertices();
  auto verts = acc_mesh->shape_vertices();
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (Int i_vert = 0; i_vert < (Int)verts.size(); ++i_vert) {
    next::Vertex::Shared_value(verts[i_vert]).set(reduction.initial_value);
  }
  #pragma omp parallel for
  for (Int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    auto& elem = elements[i_elem];
    auto& shape = elem.shape();
    for (int i_vert = 0; i_vert < nv; ++i_vert) {
      next::Vertex::Shared_value shared(shape.vertex(i_vert));
      shared.set(reduction.binary_reduction(shared.get(), get(elem, i_vert)));
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
  auto& bc_cons {_preti_masks[_preti_level]->bound_cons};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).apply_state(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_flux_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto& bc_cons {_preti_masks[_preti_level]->bound_cons};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    // write inside flux to flux cache for surface visualization/integrals
    int n_dof = params.n_dof()/params.row_size;
    Eigen::Map<Mat<>>(bc_cons[i_con].flux_cache(), n_dof) = Eigen::Map<Mat<>>(bc_cons[i_con].inside_face(true), n_dof);
    // apply boundary conditions
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).apply_flux(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_avc_diff_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto& bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).apply_diffusion(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_avc_diff_flux_bcs() {
  stopwatch["boundary conditions"].stopwatch.start();
  auto& bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).flux_diffusion(bc_cons[i_con]);
  }
  stopwatch["boundary conditions"].stopwatch.pause();
  stopwatch["boundary conditions"].work_units_completed += bc_cons.size();
}

void Solver::apply_fta_flux_bcs() {
  int rs = params.row_size;
  int nq = params.n_qpoint();
  int nv = params.n_var;
  auto& bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* in_f = bc_cons[i_con].inside_face(true);
    double* gh_f = bc_cons[i_con].ghost_face(true);
    for (int i_dof = 0; i_dof < nq*nv/rs; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
  }
}

bool Solver::use_ldg() {
  return visc.is_viscous || therm_cond.is_viscous || use_art_visc;
}

double Solver::max_dt(double msc, double msd) {
  Kernel_options opts {
    stopwatch["cartesian"],
    stopwatch["deformed"],
    stopwatch["prolong/restrict"],
    0, 0, bool(_namespace->get<int>("use_filter")),
  };
  bool local_time = _namespace->get<int>("local_time");
  if (use_ldg()) return max_dt_navier_stokes(_kernel_mesh(), opts, msc, msd, local_time, visc, therm_cond);
  else return max_dt_euler(_kernel_mesh(), opts, msc, msd, local_time);
}

void Solver::_init_face_state() {
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  auto inter {_interpreter()};
  auto& bc_cons {acc_mesh->boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh->boundary_condition(bc_sn).init_cache(bc_cons[i_con]);
    acc_mesh->boundary_condition(bc_sn).set_prescribed(inter, bc_cons[i_con]);
  }
}

Interpreter Solver::_interpreter() {
  Interpreter inter(std::vector<std::string>{});
  inter.variables = _namespace;
  return inter;
}

Solver::Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping,
               Transport_model viscosity_model, Transport_model thermal_conductivity_model,
               Turbulence_model turbulence_model,
               std::shared_ptr<Namespace> space, bool implicit)
: params{implicit ? Linearized::storage_start + Linearized::n_storage
                  : 2, n_dim + 2 + 2*(turbulence_model == k_omega),
         n_dim, row_size}
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
, _implicit{implicit}
, _preti_level{0}
{
  _namespace->assign_default("max_safety", .7); // maximum allowed safety factor for time stepping
  _namespace->assign_default("max_time_step", huge); // maximum allowed time step
  _namespace->assign_default("fix_admis_max_safety", .7); // staility ratio for fixing thermodynamic admissibility.
  _namespace->assign_default("av_diff_ratio", .3); // ratio of diffusion time to advection width
  // final scaling parameter applied to artificial viscosity coefficient
  _namespace->assign_default("av_visc_mult", 1e2);
  // maximum artificial viscosity coefficient before scaling (i.e. nondimensional)
  _namespace->assign_default("av_unscaled_max", 5.);
  _namespace->assign_default("av_advect_max_safety", .7); // stability ratio for advection
  _namespace->assign_default("av_diff_max_safety", .7); // stability ratio for diffusion
  _namespace->assign_default("buffer_dist", .8*std::sqrt(params.n_dim));
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
  _namespace->assign_default<int>("local_time", local_time_stepping);
  _namespace->assign_default("elementwise_art_visc", 0);
  _namespace->assign_default("elementwise_art_visc_diff_ratio", 5.);
  _namespace->assign_default<std::string>("working_dir", ".");
  _namespace->assign_default("iteration", 0);
  _namespace->assign_default("flow_time", 0.);
  _namespace->assign_default("time_step", 0.);
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
  calc_jacobian(false);
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

void Solver::calc_jacobian(bool snap) {
  acc_mesh->valid().assert_valid();

  // compute element jacobians
  auto& elements = acc_mesh->elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].set_jacobian(basis);
  }
  // do some extra work to make sure each face knows its normal vectors
  auto& car_cons {acc_mesh->cartesian().face_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < car_cons.size(); ++i_con) car_cons[i_con].set_normal();
  auto& def_cons {acc_mesh->deformed().face_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < def_cons.size(); ++i_con) def_cons[i_con].set_normal();
  // set position at boundary faces
  auto& bc_cons = acc_mesh->boundary_connections();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con].boundary_connection();
    HEXED_ASSERT(con.inside().element(), "connection has no element")
    Element& elem = *con.inside().element();
    con.position() = elem.face_position(basis)(con.inside().i_dim())(con.inside().sign());
    if (con.inside().is_deformed()) con.ghost().normal() = con.inside().normal();
  }
  share_vertex_data(&Element::vertex_time_step_scale, { huge, &min_fun});
  _preti_masks = acc_mesh->preti_masks(basis);
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
    for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby) {
      double s = math::chebyshev_step(n_cheby, i_cheby, cheby_safety);
      apply_avc_diff_bcs();
      opts.dt = s;
      compute_smooth_av(_kernel_mesh(), opts, [this](){apply_avc_diff_flux_bcs();}, diff_time, s);
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

  // begin estimation of high-order derivative in the style of the Cauchy-Kovalevskaya theorem using a linear advection equation.
  // perform pseudotime iteration
  for (int iter = 0; iter < _namespace->get<int>("av_advect_iters"); ++iter) {
    sw_adv["setup"].stopwatch.start();
    // evaluate advection operator
    compute_write_face_advection(_kernel_mesh());
    compute_prolong_advection(_kernel_mesh());
    sw_adv["setup"].stopwatch.pause();
    for (int i = 0; i < 2; ++i) {
      sw_adv["BCs"].stopwatch.start();
      auto& bc_cons {acc_mesh->boundary_connections()};
      #pragma omp parallel for
      for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
        int bc_sn = bc_cons[i_con].bound_cond_serial_n();
        acc_mesh->boundary_condition(bc_sn).apply_advection(bc_cons[i_con]);
      }
      sw_adv["BCs"].stopwatch.pause();
      sw_adv["BCs"].work_units_completed += acc_mesh->elements().size();
      opts.i_stage = i;
      compute_advection(_kernel_mesh(), opts, advect_length);
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
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double proj = 0;
      for (int i_proj = 0; i_proj < rs; ++i_proj) {
        proj += adv[i_proj*nq + i_qpoint]*weights(i_proj)*orth(i_proj);
      }
      forcing[i_qpoint] = proj*proj*2*state[(nd + 1)*nq + i_qpoint]/state[nd*nq + i_qpoint];
    }
  } // Cauchy-Kovalevskaya-style derivative estimate complete!

  // begin root-smear-square operation
  int n_real = params.n_forcing - 1; // number of real time steps (as apposed to pseudotime steps)
  double diff_time = _namespace->get<double>("av_diff_ratio")*advect_length*advect_length/n_real; // compute size of real time step (as opposed to pseudotime)
  stopwatch["set art visc"]["diffusion"].stopwatch.start();
  diffuse_art_visc(diff_time);
  stopwatch["set art visc"]["diffusion"].stopwatch.pause();
  stopwatch["set art visc"]["diffusion"].work_units_completed += elements.size();

  // clean up
  double mult = _namespace->get<double>("av_visc_mult")*advect_length;
  double us_max = advect_length*_namespace->get<double>("av_unscaled_max")*std::sqrt(2*_namespace->get<double>("freestream" + std::to_string(nd + 1))/_namespace->get<double>("freestream" + std::to_string(nd)));
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
                      {-huge, &max_fun});
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
                    {-huge, &max_fun});
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
    printers::info(to_string(elem.uncertainty) + "\n");
  }
}

void Solver::update() {
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  double safety = _namespace->get<double>("max_safety");
  double cheby_safety = _namespace->get<double>("cheby_safety");
  for (int i_flow = 0; i_flow < _namespace->get<int>("flow_iters"); ++i_flow) {
    // compute time step
    double dt = 0;
    HEXED_ASSERT(_preti_masks.size(), "meshing mask list is empty");
    int n_preti = (_namespace->get<int>("bl_multirate") && !i_flow) ? _preti_masks.size() : 1;
    for (int i_preti = 0; i_preti < n_preti; ++i_preti) if (i_preti != 1) {
      int n_bl = i_preti ? _namespace->get<int>("bl_iters") : 1;
      for (int i_bl = 0; i_bl < n_bl; ++i_bl) {
        int n_cheby = i_preti ? _namespace->get<int>("n_cheby_bl") : _namespace->get<int>("n_cheby_flow");
        int max_sub_iters = i_preti ? _namespace->get<int>("max_conv_sub_iters") : 1;
        double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1, cheby_safety);
        // run chebyshev iterations
        for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby) {
          _preti_level = i_preti - bool(i_preti);
          Kernel_mesh& km = _preti_masks[_preti_level]->kernel_mesh;
          double cheby_step = math::chebyshev_step(n_cheby, i_cheby, cheby_safety);
          int sub_iters = std::ceil(max_sub_iters*cheby_step/max_cheby - 1e-6);
          double nominal_dt = std::min(max_dt(safety/max_cheby*sub_iters, safety), _namespace->get<double>("max_time_step"));
          dt = nominal_dt*cheby_step;
          HEXED_ASSERT(!std::isnan(dt), "time step is NaN", assert::Numerical_exception);
          bool fixed = false;
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
              };
              apply_state_bcs();
              if (use_ldg() && !i && !i_sub) compute_navier_stokes(km, opts, [this](){apply_flux_bcs();}, visc, therm_cond, _namespace->get<int>("iteration")%100000 == 0 && _namespace->get<int>("iteration") != 0);
              else compute_euler(km, opts);
              // note that function call must come first to ensure it is evaluated despite short-circuiting
              fixed = fix_admissibility(_namespace->get<double>("fix_admis_max_safety")) || fixed;
            }
            stopwatch.work_units_completed += km.elems.size();
            stopwatch["cartesian"].work_units_completed += km.car_elems.size();
            stopwatch["deformed" ].work_units_completed += km.def_elems.size();
            if (fixed) break;
          }
          if (fixed) break;

          // update status for reporting
          _namespace->assign<double>("time_step", dt);
          _namespace->assign<double>("flow_time", _namespace->get<double>("flow_time") + dt);
          status.time_step = dt;
          status.flow_time += dt;
        }
      }
    }
  }
  _preti_level = 0;

  ++status.iteration;
  stopwatch.stopwatch.pause();
}

void Solver::update_implicit() {
  HEXED_ASSERT(_implicit, "`update_implicit` called on a Solver that was not constructed in implicit mode");
  Linearized lin(*this);
  iterative::gmres(lin, 27, 1);
  lin.add(-Linearized::storage_start, 1., -Linearized::storage_start + 3, 1., 0.);
  compute_write_face(_kernel_mesh());
  compute_prolong(_kernel_mesh());
  fix_admissibility(.7);
}

void Solver::compute_residual() {
  apply_state_bcs();
  Kernel_options opts {
    stopwatch["cartesian"],
    stopwatch["deformed"],
    stopwatch["prolong/restrict"],
    1.,
    0,
    true,
    bool(_namespace->get<int>("use_filter")),
  };
  if (use_ldg()) compute_navier_stokes(_kernel_mesh(), opts, [this](){apply_flux_bcs();}, visc, therm_cond, false);
  else compute_euler(_kernel_mesh(), opts);
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
    max_dt(safeties[i_term]/max_cheby, safeties[!i_term]);
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
  bool admiss = 1;
  auto check_admis = [&](double* data, int n_qpoint, int n_var) {
    bool adm = true;
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      adm = adm && (data[nd*n_qpoint + i_qpoint] > 0.)
                && (data[(nd + 1)*n_qpoint + i_qpoint] > 0.);
      for (int i_var = 0; i_var < n_var; ++i_var) {
        HEXED_ASSERT(std::isfinite(data[i_var*n_qpoint + i_qpoint]),
                     format_str(200, "variable %i = %e has non-finite value.", i_var, data[i_var*n_qpoint + i_qpoint]),
                     assert::Numerical_exception);
      }
    }
    return adm;
  };
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].record = 0;
  }
  #pragma omp parallel for reduction (&&:admiss)
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
  auto& ref_faces = _preti_masks[_preti_level]->kernel_mesh.ref_faces;
  bool refined_admiss = 1;
  #pragma omp parallel for reduction (&&:refined_admiss)
  for (int i_face = 0; i_face < ref_faces.size(); ++i_face) {
    auto& ref = ref_faces[i_face];
    int n_fine = params.n_vertices()/2;
    for (int i_dim = 0; i_dim < nd - 1; ++i_dim) n_fine /= 1 + ref.stretch[i_dim];
    for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
      refined_admiss = refined_admiss && check_admis(ref.fine[i_fine], nq/rs, nd + 2);
    }
  }
  sw.work_units_completed += acc_mesh->elements().size();
  sw.stopwatch.pause();
  return admiss && refined_admiss;
}

bool Solver::fix_admissibility(double stability_ratio) {
  if (!fix_admis) return false;
  auto& sw_fix = stopwatch["fix admis."];
  sw_fix.stopwatch.start();
  std::string wd = _namespace->get<std::string>("working_dir");
  std::string vis_expr = _namespace->get<std::string>("vis_field_vars");
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  const int nv = params.n_vertices();
  int iter;
  int n_iters = std::numeric_limits<int>::max();
  for (iter = 0; iter < n_iters; ++iter) {
    HEXED_ASSERT(iter < 1e5, format_str(200, "failed to fix thermodynamic admissability in %i iterations", iter));
    if (is_admissible()) {
      if (iter) n_iters = std::min(n_iters, 2*iter);
      else {
        ++iter;
        break;
      }
    } else {
      n_iters = std::numeric_limits<int>::max();
    }
    if (iter == 100) visualize_field("default", wd + "severe_indamis" + std::to_string(status.iteration), vis_expr);
    if (iter == 0) {
      printers::warn("Warning: ", true);
      printers::warn(format_str(200, "Thermodynamically inadmissible state detected (solver iteration %i). Attempting to fix...\n",
                                _namespace->get<int>("iteration")));
    }
    printers::warn(format_str(200, "    iteration %i\n", iter));
    auto& elems = acc_mesh->elements();
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        elem.vertex_fix_admis_coef(i_vert) = elem.record;
      }
    }
    share_vertex_data(&Element::vertex_fix_admis_coef, {-huge, &max_fun});
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      double max_fac = 0;
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        max_fac = std::max(max_fac, elem.vertex_fix_admis_coef(i_vert));
      }
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        elem.vertex_fix_admis_coef(i_vert) = max_fac;
      }
    }
    share_vertex_data(&Element::vertex_fix_admis_coef, {-huge, &max_fun});
    Mat<dyn, dyn> interp(rs, 2);
    interp(all, 0) = Mat<>::Ones(rs) - basis.nodes();
    interp(all, 1) = basis.nodes();
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      Mat<> vert_fac(nv);
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        vert_fac(i_vert) = elem.vertex_fix_admis_coef(i_vert);
      }
      Eigen::Map<Mat<>>(elem.laplacian_av_coef(), nq) = math::hypercube_matvec(interp, vert_fac);
    }
    if (status.iteration >= last_fix_vis_iter + 1000 && iter == 0) {
      last_fix_vis_iter = status.iteration;
      visualize_field("default", wd + "inadmis" + std::to_string(status.iteration), vis_expr);
    }
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        std::swap(elem.laplacian_av_coef()[i_qpoint], elem.bulk_av_coef()[i_qpoint]);
      }
    }
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
    dt = 1.;
    double linear = dt;
    double quadratic = dt*dt/8/0.9;
    std::array<double, 2> step;
    step[1] = (linear + std::sqrt(linear*linear - 4*quadratic))/2.;
    step[0] = quadratic/step[1];
    for (double s : step) {
      auto& bc_cons {acc_mesh->boundary_connections()};
      #pragma omp parallel for
      for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
        double* in_f = bc_cons[i_con].inside_face(false);
        double* gh_f = bc_cons[i_con].ghost_face(false);
        for (int i_dof = 0; i_dof < nq*params.n_var/rs; ++i_dof) gh_f[i_dof] = in_f[i_dof];
      }
      opts.dt = s;
      compute_fix_therm_admis(_kernel_mesh(), opts, [this](){apply_fta_flux_bcs();});
    }
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        std::swap(elem.laplacian_av_coef()[i_qpoint], elem.bulk_av_coef()[i_qpoint]);
      }
    }
  }
  --iter;
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

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func& func) {
  return func(acc_mesh->element(ref_level, is_deformed, serial_n), basis, i_qpoint, _namespace->get<double>("flow_time"));
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
    double volume = math::pow(element.nominal_size(), params.n_dim);
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
  Vis_evaluator(Interpreter&& inter, std::function<void(Namespace&, T&)> assign, std::string expr, T& t, int n_dim_topo)
  : _inter{inter}, _assign{assign}, _expr{expr}, _n_dim_topo{n_dim_topo}
  {
    Storage_params params = t.storage_params();
    _n_dim = params.n_dim;
    auto sub = _inter.make_sub();
    _assign(*sub.variables, t);
    sub.subspace();
    sub.exec(_expr);
    _var_names = sub.variables->names();
    _n_var = _var_names.size();
    _shape = hypercubes(_n_dim + _n_var, _n_dim_topo, params.row_size);
  }

  Array<double> evaluate(T& t) {
    Array<double> qpoints(_shape);
    auto sub = _inter.make_sub();
    _assign(*sub.variables, t);
    sub.exec(_expr);
    for (int i_dim = 0; i_dim < _n_dim; ++i_dim) {
      qpoints(i_dim) = sub.variables->get<Array<double>>("pos" + std::to_string(i_dim));
    }
    for (int i_var = 0; i_var < _n_var; ++i_var) {
      sub.variables->assign_array(qpoints(_n_dim + i_var), _var_names[i_var]);
    }
    return qpoints;
  }
  std::vector<std::string> var_names() {return _var_names;}

  void visualize(std::string format, std::string name, int n_sample, bool wireframe, Sequence<T&>& seq,
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
  auto& bc_cons {acc_mesh->boundary_connections()};
  if (!bc_cons.size()) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, bc_cons[0], params.n_dim - 1
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
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con {bc_cons[i_con]};
    if (con.bound_cond_serial_n() != bc_sn) continue;
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
  auto& bc_cons {acc_mesh->boundary_connections()};
  if (!bc_cons.size()) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, bc_cons[0], params.n_dim - 1
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
    if (con.bound_cond_serial_n() != bc_sn) continue;
    auto& elem = con.element();
    double area = math::pow(elem.nominal_size(), nd - 1);
    Array<double> qpoints{evaluator.evaluate(con)};
    Kernel_connection& ker_con = con;
    double* nrml = ker_con.normal();
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      double nrml_mag = 0;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        nrml_mag += math::pow(nrml[i_dim*nfq + i_qpoint], 2);
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
  evaluator.visualize(format, name, n_sample, wireframe, elems,
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
  auto& bc_cons {acc_mesh->boundary_connections()};
  if (!bc_cons.size()) return;
  Vis_evaluator<Boundary_connection> evaluator(
    _interpreter(),
    [&](Namespace& space, Boundary_connection& con){vis_variables::surface(space, con);},
    expr, bc_cons[0], params.n_dim - 1
  );
  evaluator.visualize(format, name, n_sample, wireframe, bc_cons,
                      _namespace->get<double>("flow_time"), basis,
                      [bc_sn](Boundary_connection& con){return con.bound_cond_serial_n() == bc_sn;});
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
    max_dt(safeties[i_term], safeties[!i_term]);
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
