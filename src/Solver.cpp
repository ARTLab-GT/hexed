#include <fstream>
#include <iostream>

#include <config.hpp>
#include <Solver.hpp>
#include <Tecplot_file.hpp>
#include <Vis_data.hpp>
#include <thermo.hpp>
#include <Xdmf_wrapper.hpp>
#include <iterative.hpp>
#include <Gauss_lobatto.hpp>
#include <kernels.hpp>

// kernels
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>
#include <Face_permutation.hpp>

#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed
{

void Solver::share_vertex_data(std::function<double&(Element&, int i_vertex)> access_fun, std::function<double(Mat<>)> reduce)
{
  share_vertex_data(access_fun, access_fun, reduce);
}

void Solver::share_vertex_data(std::function<double(Element&, int i_vertex)> get, std::function<double&(Element&, int i_vertex)> set, std::function<double(Mat<>)> reduce)
{
  auto& elements = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].push_shareable_value(get);
  }
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].fetch_shareable_value(set, reduce);
  }
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  #pragma omp parallel for
  for (int i_match = 0; i_match < matchers.size(); ++i_match) matchers[i_match].match(set);
}

void Solver::apply_state_bcs()
{
  stopwatch.children.at("boundary conditions").stopwatch.start();
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->apply_state(bc_cons[i_con]);
  }
  stopwatch.children.at("boundary conditions").stopwatch.pause();
  stopwatch.children.at("boundary conditions").work_units_completed += bc_cons.size();
}

void Solver::apply_flux_bcs()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    // write inside flux to flux cache for surface visualization/integrals
    int n_dof = params.n_dof()/params.row_size;
    Eigen::Map<Mat<>>(bc_cons[i_con].flux_cache(), n_dof) = Eigen::Map<Mat<>>(bc_cons[i_con].inside_face() + 2*n_dof, n_dof);
    // apply boundary conditions
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->apply_flux(bc_cons[i_con]);
  }
}

void Solver::apply_avc_diff_bcs()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->apply_diffusion(bc_cons[i_con]);
  }
}

void Solver::apply_avc_diff_flux_bcs()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->flux_diffusion(bc_cons[i_con]);
  }
}

void Solver::apply_fta_flux_bcs()
{
  int nd = params.n_dim;
  int rs = params.row_size;
  int nq = params.n_qpoint();
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* in_f = bc_cons[i_con].inside_face() + 2*(nd + 2)*nq/rs;
    double* gh_f = bc_cons[i_con].ghost_face()  + 2*(nd + 2)*nq/rs;
    for (int i_dof = 0; i_dof < nq*(nd + 2)/rs; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
  }
}

bool Solver::use_ldg()
{
  return visc.is_viscous || therm_cond.is_viscous || use_art_visc;
}

double Solver::max_dt(double msc, double msd)
{
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  double use_filt = _namespace->lookup<int>("use_filter").value();
  bool local_time = _namespace->lookup<int>("local_time").value();
  #define MAX_DT(is_visc, ...) { \
    return std::min((*kernel_factory<Spatial<pde::Navier_stokes<is_visc>::Pde, false>::Max_dt>(nd, rs, basis, local_time, use_filt, msc, msd, __VA_ARGS__))(acc_mesh.cartesian().kernel_elements(), sw_car, "compute time step"), \
                    (*kernel_factory<Spatial<pde::Navier_stokes<is_visc>::Pde,  true>::Max_dt>(nd, rs, basis, local_time, use_filt, msc, msd, __VA_ARGS__))(acc_mesh.deformed ().kernel_elements(), sw_def, "compute time step")); \
  }
  if (use_ldg()) MAX_DT(true , visc, therm_cond)
  else           MAX_DT(false, visc, therm_cond)
  #undef MAX_DT
}

Solver::Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping,
               Transport_model viscosity_model, Transport_model thermal_conductivity_model,
               std::shared_ptr<Namespace> space, std::shared_ptr<Printer> printer, bool implicit) :
  params{implicit ? Linearized::storage_start + Linearized::n_storage : 2, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size},
  basis{row_size},
  stopwatch{"(element*iteration)"},
  use_art_visc{false},
  fix_admis{false},
  av_rs{row_size},
  write_face(kernel_factory<Spatial<pde::Navier_stokes<false>::Pde, false>::Write_face>(params.n_dim, params.row_size, basis)), // note: for now, false and true are equivalent for `Write_face`
  visc{viscosity_model},
  therm_cond{thermal_conductivity_model},
  _namespace{space},
  _printer{printer},
  _implicit{implicit}
{
  _namespace->assign_default("max_safety", .7); // maximum allowed safety factor for time stepping
  _namespace->assign_default("max_time_step", huge); // maximum allowed time step
  _namespace->assign_default("fix_admis_max_safety", .7); // staility ratio for fixing thermodynamic admissibility.
  _namespace->assign_default("av_diff_ratio", .3); // ratio of diffusion time to advection width
  _namespace->assign_default("av_visc_mult", 1e2); // final scaling parameter applied to artificial viscosity coefficient
  _namespace->assign_default("av_unscaled_max", 5.); // maximum artificial viscosity coefficient before scaling (i.e. nondimensional)
  _namespace->assign_default("av_advect_max_safety", .7); // stability ratio for advection
  _namespace->assign_default("av_diff_max_safety", .7); // stability ratio for diffusion
  _namespace->assign_default("buffer_dist", .8*std::sqrt(params.n_dim));
  _namespace->assign_default("n_cheby_flow", 1);
  _namespace->assign_default("n_cheby_av", 1);
  _namespace->assign_default("av_advect_iters", 1); // number of advection iterations to run each time `update_art_visc_smoothness` is called
  _namespace->assign_default("av_diff_iters", 1); // number of diffusion iterations to run each time `update_art_visc_smoothness` is called
  _namespace->assign_default("flow_iters", 1);
  _namespace->assign_default("use_filter", 0); // whether to use modal filter acceleration
  _namespace->assign_default<int>("local_time", local_time_stepping);
  _namespace->assign_default("elementwise_art_visc", 0);
  _namespace->assign_default("elementwise_art_visc_diff_ratio", 5.);
  _namespace->assign_default("laplacian_art_visc", 0);
  _namespace->assign_default<std::string>("working_dir", ".");
  _namespace->assign("fix_iters", 0);
  _namespace->assign("iteration", 0);
  _namespace->assign("flow_time", 0.);
  _namespace->assign("av_advection_residual", 0.);
  _namespace->assign("av_diffusion_residual", 0.);
  _namespace->assign("wall_time", 0.);
  status.set_time();
  // setup categories for performance reporting
  std::string unit = "(element*(time integration stage))";
  stopwatch.children.emplace("prolong/restrict", unit);
  stopwatch.children.emplace("fix admis.", "(element*(fix admis. iter))");
  stopwatch.children.emplace("check admis.", "(element*iteration)");
  stopwatch.children.emplace("set art visc", stopwatch.work_unit_name);
  stopwatch.children.at("set art visc").children.emplace("initialize", stopwatch.work_unit_name);
  stopwatch.children.at("set art visc").children.emplace("advection", stopwatch.work_unit_name);
  stopwatch.children.at("set art visc").children.at("advection").children.emplace("update", unit);
  stopwatch.children.at("set art visc").children.at("advection").children.emplace("setup", stopwatch.work_unit_name);
  stopwatch.children.at("set art visc").children.at("advection").children.emplace("BCs", unit);
  stopwatch.children.at("set art visc").children.emplace("diffusion", stopwatch.work_unit_name);
  for (std::string type : {"cartesian", "deformed"}) {
    for (auto* sw : {&stopwatch, &stopwatch.children.at("set art visc").children.at("advection"),
                     &stopwatch.children.at("set art visc").children.at("diffusion"), &stopwatch.children.at("fix admis.")}) {
      sw->children.emplace(type, stopwatch.work_unit_name);
      auto& children = sw->children.at(type).children;
      children.emplace("compute time step", stopwatch.work_unit_name);
      children.emplace("neighbor", "(connection*(time integration stage))");
      children.emplace("local", unit);
    }
    for (auto* sw : {&stopwatch, &stopwatch.children.at("fix admis."), &stopwatch.children.at("set art visc").children.at("diffusion")}) {
      sw->children.at(type).children.emplace("reconcile LDG flux", unit);
    }
  }
  stopwatch.children.emplace("boundary conditions", "(boundary connection)*(time integration stage)");
  // initialize advection state to 1
  auto& elements = acc_mesh.elements();
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

Mesh& Solver::mesh() {return acc_mesh;}
Storage_params Solver::storage_params() {return params;}
const Stopwatch_tree& Solver::stopwatch_tree() {return stopwatch;}

void Solver::snap_faces()
{
  const int nd = params.n_dim;
  const int nfq = params.n_qpoint()/params.row_size;
  auto& bc_cons {acc_mesh.boundary_connections()};
  auto& elems = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Lock::Acquire acq(bc_cons[i_con].element().lock);
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).mesh_bc->snap_node_adj(bc_cons[i_con], basis);
  }
  Gauss_lobatto lob(std::max(2, basis.row_size - 1));
  Mat<dyn, dyn> to_lob = basis.interpolate(lob.nodes());
  Mat<dyn, dyn> from_lob = lob.interpolate(basis.nodes());
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].stage(0);
    double* rk_ref = elems[i_elem].stage(1);
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) rk_ref[i_dof] = state[i_dof];
  }
  for (int sweep = 0; sweep < nd - 1; ++sweep) {
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      double* state = elems[i_elem].stage(0);
      for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = 0;
    }
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      int bc_sn = con.bound_cond_serial_n();
      if (bc_sn == acc_mesh.surface_bc_sn()) {
        double* state = con.element().stage(0);
        double* node_adj = con.element().node_adjustments();
        if (node_adj) {
          node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
          for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
            for (int i_node = 0; i_node < params.row_size; ++i_node) {
              state[ind.i_qpoint(i_node)] = math::sign(con.inside_face_sign())*node_adj[ind.i_face_qpoint()];
            }
          }
        }
      }
    }
    (*write_face)(acc_mesh.kernel_elements());
    (*kernel_factory<Prolong_refined>(params.n_dim, params.row_size, basis))(acc_mesh.refined_faces());
    Copy fake_bc;
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      fake_bc.apply_state(con);
    }
    if (sweep) {
      auto& ref_cons = acc_mesh.deformed().refined_connections();
      #pragma omp parallel for
      for (int i_con = 0; i_con < ref_cons.size(); ++i_con) {
        for (int i_fine = 0; i_fine < 2; ++i_fine) {
          double* state = ref_cons[i_con].connection(i_fine).state();
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
            state[i_qpoint] *= 2;
            state[nfq*params.n_var + i_qpoint] = 0;
          }
        }
      }
    }
    (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Neighbor>(params.n_dim, params.row_size, 0))(acc_mesh.deformed().kernel_connections());
    (*kernel_factory<Restrict_refined>(params.n_dim, params.row_size, basis, false, true))(acc_mesh.refined_faces());
    #pragma omp parallel for
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      auto& con = bc_cons[i_con];
      Lock::Acquire acq(con.element().lock);
      int bc_sn = con.bound_cond_serial_n();
      if (bc_sn == acc_mesh.surface_bc_sn()) {
        double* state = con.element().stage(0);
        double* node_adj = con.element().node_adjustments();
        if (node_adj) {
          node_adj += (2*con.i_dim() + con.inside_face_sign())*params.n_qpoint()/params.row_size;
          for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) if (j_dim != con.i_dim()) {
            for (Row_index ind(params.n_dim, params.row_size, j_dim); ind; ++ind) {
              Mat<> row(params.row_size);
              for (int i_node = 0; i_node < params.row_size; ++i_node) {
                row(i_node) = state[ind.i_qpoint(i_node)];
              }
              Mat<> lrow = to_lob*row;
              for (int i_sign = 0; i_sign < 2; ++i_sign) {
                lrow(i_sign*(lob.row_size - 1)) = con.element().faces[2*j_dim + i_sign][2*params.n_var*params.n_qpoint()/params.row_size + ind.i_face_qpoint()];
              }
              row = from_lob*lrow;
              for (int i_node = 0; i_node < params.row_size; ++i_node) {
                state[ind.i_qpoint(i_node)] = row(i_node);
              }
            }
          }
          for (Row_index ind(params.n_dim, params.row_size, con.i_dim()); ind; ++ind) {
            node_adj[ind.i_face_qpoint()] = math::sign(con.inside_face_sign())*state[ind.i_qpoint(0)];
          }
        }
      }
    }
  }
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].stage(0);
    double* rk_ref = elems[i_elem].stage(1);
    for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) state[i_dof] = rk_ref[i_dof];
  }
}

void Solver::calc_jacobian()
{
  acc_mesh.valid().assert_valid();
  snap_faces();
  const int n_dim = params.n_dim;
  const int rs = basis.row_size;
  const int nfq = params.n_qpoint()/rs;
  const int nfdof = nfq*params.n_var;

  // compute element jacobians
  auto& elements = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].set_jacobian(basis);
  }

  /*
   * compute surface normals for deformed connections
   */
  auto& def_cons {acc_mesh.deformed().face_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < def_cons.size(); ++i_con) {
    double* nrml = def_cons[i_con].normal();
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) nrml[i_data] = 0.;
  }
  // for deformed refined faces, set normal to coarse face normal (for Cartesian, setting normal is not necessary)
  auto& ref_cons = acc_mesh.deformed().refined_connections();
  (*kernel_factory<Prolong_refined>(n_dim, rs, basis, true))(acc_mesh.refined_faces());
  #pragma omp parallel for
  for (int i_ref = 0; i_ref < ref_cons.size(); ++i_ref) {
    auto& ref = ref_cons[i_ref];
    bool rev = ref.order_reversed();
    auto dir = ref.direction();
    // might need to flip the direction of the normal vector depending on connection direction
    int sign = 1 - 2*(dir.flip_normal(0) != dir.flip_normal(1));
    for (int i_fine = 0; i_fine < ref.n_fine_elements(); ++i_fine) {
      auto& fine = ref.connection(i_fine);
      double* face [2] {fine.state() + rev*nfdof, fine.state() + (!rev)*nfdof};
      auto fp = kernel_factory<Face_permutation>(n_dim, rs, dir, face[1]);
      fp->match_faces();
      for (int i_data = 0; i_data < n_dim*nfq; ++i_data) {
        face[1][i_data] = sign*face[0][i_data];
      }
      fp->restore();
    }
  }
  // for BCs, copy normal to ghost face
  auto& bc_cons = acc_mesh.boundary_connections();
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* in_f = bc_cons[i_con].inside_face();
    double* gh_f = bc_cons[i_con].ghost_face();
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) gh_f[i_data] = in_f[i_data];
  }
  // compute the shared face normal
  #pragma omp parallel for
  for (int i_con = 0; i_con < def_cons.size(); ++i_con)
  {
    auto& con = def_cons[i_con];
    double* elem_nrml [2] {con.state(), con.state() + nfdof};
    auto dir = con.direction();
    // permute face 1 so that quadrature points match up
    auto fp = kernel_factory<Face_permutation>(n_dim, rs, dir, elem_nrml[1]);
    fp->match_faces();
    // take average of element face normals with appropriate flipping
    int sign [2];
    for (int i_side : {0, 1}) sign[i_side] = 1 - 2*dir.flip_normal(i_side);
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) {
      double n = 0;
      for (int i_side : {0, 1}) n += 0.5*sign[i_side]*elem_nrml[i_side][i_data];
      for (int i_side : {0, 1}) elem_nrml[i_side][i_data] = sign[i_side]*n;
    }
    // put face 1 back in its original order (except now we're working with the normal data not state data)
    fp->restore();
    for (int i_side = 0; i_side < 2; ++i_side) {
      double* n = con.normal(i_side);
      for (int i_data = 0; i_data < n_dim*nfq; ++i_data) {
        n[i_data] = elem_nrml[i_side][i_data];
      }
    }
  }
  // write face normal for coarse hanging node faces
  #pragma omp parallel for
  for (int i_ref = 0; i_ref < ref_cons.size(); ++i_ref) {
    auto& ref = ref_cons[i_ref];
    bool rev = ref.order_reversed();
    auto& elem = ref.connection(0).element(rev);
    auto dir = ref.direction();
    int i_face = 2*dir.i_dim[rev] + dir.face_sign[rev];
    double* nrml = elem.face_normal(i_face);
    double* state = elem.faces[i_face];
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) {
      nrml[i_data] = state[i_data];
    }
  }
  // set position at boundary faces
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Element& elem = bc_cons[i_con].element();
    double* surf_pos = bc_cons[i_con].surface_position();
    int i_face = bc_cons[i_con].direction().i_face(0);
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      auto pos = elem.face_position(basis, i_face, i_qpoint);
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        surf_pos[i_dim*nfq + i_qpoint] = pos[i_dim];
      }
    }
  }
  share_vertex_data(&Element::vertex_time_step_scale, Vertex::vector_min);
}

void Solver::initialize(const Spacetime_func& func)
{
  acc_mesh.valid().assert_valid();
  if (func.n_var(params.n_dim) < params.n_var) {
    throw std::runtime_error("initializer has too few output variables");
  }
  auto& elements = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
      std::vector<double> pos_vec {};
      auto state = func(elements[i_elem].position(basis, i_qpoint), status.flow_time);
      for (int i_var = 0; i_var < params.n_var; ++i_var) {
        elements[i_elem].stage(0)[i_var*params.n_qpoint() + i_qpoint] = state[i_var];
      }
    }
  }
  (*write_face)(acc_mesh.kernel_elements());
  (*kernel_factory<Prolong_refined>(params.n_dim, params.row_size, basis))(acc_mesh.refined_faces());
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->init_cache(bc_cons[i_con]);
  }
}

void Solver::set_art_visc_off()
{
  use_art_visc = false;
}

void Solver::set_art_visc_constant(double value)
{
  use_art_visc = true;
  auto& elements = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* av = elements[i_elem].art_visc_coef();
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
     av[i_qpoint] = value;
    }
  }
}

void Solver::diffuse_art_visc(int n_real, double diff_time)
{
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  auto& elements = acc_mesh.elements();
  // evaluate CFL condition
  double diff_safety = _namespace->lookup<double>("av_diff_max_safety").value();
  double n_cheby = _namespace->lookup<double>("n_cheby_av").value();
  (*kernel_factory<Spatial<pde::Smooth_art_visc, false>::Max_dt>(nd, rs, basis, true, false, 1., diff_safety))(acc_mesh.cartesian().kernel_elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("cartesian"), "compute time step");
  (*kernel_factory<Spatial<pde::Smooth_art_visc,  true>::Max_dt>(nd, rs, basis, true, false, 1., diff_safety))(acc_mesh.deformed ().kernel_elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("deformed" ), "compute time step");
  // initialize residual to zero (will compute RMS over all real time steps)
  status.diff_res = 0;
  for (int real_step = 0; real_step < n_real; ++real_step)
  {
    double diff = 0;
    int n_avg = 0;
    // copy value from forcing function storage to scalar variable of stage 0
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      double* state = elements[i_elem].stage(0);
      double* forcing = elements[i_elem].art_visc_forcing();
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        state[i_qpoint] = forcing[(real_step + 1)*nq + i_qpoint];
      }
    }
    (*write_face)(acc_mesh.kernel_elements());
    (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
    diff = 0;
    n_avg = 0;
    // perform pseudotime iteration
    for (int i_iter = 0; i_iter < _namespace->lookup<int>("av_diff_iters").value(); ++i_iter)
    {
      // record initial state for residual calculation
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
        double* state = elements[i_elem].stage(0);
        double* forcing = elements[i_elem].art_visc_forcing();
        for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
          forcing[(real_step + 1)*nq + i_qpoint] = state[i_qpoint];
        }
      }
      for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby) {
        double s = math::chebyshev_step(n_cheby, i_cheby);
        apply_avc_diff_bcs();
        // update state with spatial term
        compute_avc_diff(s, 0);
        // update state with real time forcing term
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* state = elements[i_elem].stage(0);
          double* forcing = elements[i_elem].art_visc_forcing();
          double* tss = elements[i_elem].time_step_scale();
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
            double pseudotime_scale = s/diff_time*tss[i_qpoint];
            double f = forcing[real_step*nq + i_qpoint];
            f = (real_step == 1) ? std::sqrt(f) : f; // at time step 1 switch from square root domain to linear domain
            state[i_qpoint] += f*pseudotime_scale;
            state[i_qpoint] /= 1 + pseudotime_scale;
          }
        }
        // update face state
        (*write_face)(acc_mesh.kernel_elements());
        (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
      }
    }
    // update forcing function and residual
    #pragma omp parallel for reduction(+:diff, n_avg)
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      double* state = elements[i_elem].stage(0);
      double* forcing = elements[i_elem].art_visc_forcing();
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        // update residual
        double d = forcing[(real_step + 1)*nq + i_qpoint] - state[i_qpoint];
        diff += d*d;
        ++n_avg;
        // update forcing
        forcing[(real_step + 1)*nq + i_qpoint] = std::max(0., state[i_qpoint]);
      }
    }
    status.diff_res += diff/n_avg;
  }
  status.diff_res = std::sqrt(status.diff_res/n_real); // finish computing RMS residual
  _namespace->assign("av_diffusion_residual", std::sqrt(status.diff_res/n_real));
}

void Solver::update_art_visc_smoothness(double advect_length)
{
  stopwatch.stopwatch.start();
  stopwatch.children.at("set art visc").stopwatch.start();
  use_art_visc = true;
  const int nq = params.n_qpoint();
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& elements = acc_mesh.elements();

  // store the current state in stage 1 (normally the time integration reference)
  // so we can use stage 0 for solving the advection equation
  stopwatch.children.at("set art visc").children.at("initialize").stopwatch.start();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].stage(0);
    double* rk_ref = elements[i_elem].stage(1);
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      for (int i_var = 0; i_var < params.n_var; ++i_var) {
        int i = i_var*nq + i_qpoint;
        rk_ref[i] = state[i];
      }
    }
  }
  // set advection velocity
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].stage(0);
    double* rk_ref = elements[i_elem].stage(1);
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double scale = sqrt(2*rk_ref[nd*nq + i_qpoint]*rk_ref[(nd + 1)*nq + i_qpoint]);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        state[i_dim*nq + i_qpoint] = rk_ref[i_dim*nq + i_qpoint]/scale;
      }
    }
  }
  stopwatch.children.at("set art visc").children.at("initialize").stopwatch.pause();
  stopwatch.children.at("set art visc").children.at("initialize").work_units_completed += elements.size();
  // enforce CFL condition
  auto& sw_adv = stopwatch.children.at("set art visc").children.at("advection");
  sw_adv.stopwatch.start();
  (*write_face)(acc_mesh.kernel_elements());
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
  double adv_safety = _namespace->lookup<double>("av_advect_max_safety").value();
  (*kernel_factory<Spatial<pde::Advection, false>::Max_dt>(nd, rs, basis, true, false, adv_safety, 1.))(acc_mesh.cartesian().kernel_elements(), sw_adv.children.at("cartesian"), "compute time step");
  (*kernel_factory<Spatial<pde::Advection,  true>::Max_dt>(nd, rs, basis, true, false, adv_safety, 1.))(acc_mesh.deformed ().kernel_elements(), sw_adv.children.at("deformed" ), "compute time step");

  // begin estimation of high-order derivative in the style of the Cauchy-Kovalevskaya theorem using a linear advection equation.
  double diff = 0; // for residual computation
  int n_avg = 0;
  // perform pseudotime iteration
  for (int iter = 0; iter < _namespace->lookup<int>("av_advect_iters").value(); ++iter)
  {
    diff = 0;
    n_avg = 0;
    for (int sign : {-1, 1}) // do forward and backward advection separately
    {
      // flip direction of advection velocity
      #pragma omp parallel for
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
        double* state = elements[i_elem].stage(0);
        for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
          for (int i_dim = 0; i_dim < nd; ++i_dim) state[i_dim*nq + i_qpoint] *= -1;
        }
      }
      // for odd row size we have to set the advection state corresponding to t~ = 0 to 1
      if (rs%2) {
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* adv = elements[i_elem].advection_state();
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) adv[(rs/2)*nq + i_qpoint] = 1;
        }
      }
      // loop through nodes of projection basis
      for (int i_proj = rs - rs/2; i_proj < rs; ++i_proj)
      {
        // find out which node we're advecting to
        int i_node = (sign > 0) ? i_proj : rs - 1 - i_proj; // if sign < 0 we are looping through the nodes backward
        // copy advection state to the scalar variables of stage 1
        sw_adv.children.at("setup").stopwatch.start();
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* state = elements[i_elem].stage(0);
          double* adv = elements[i_elem].advection_state();
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) state[nd*nq + i_qpoint] = adv[i_node*nq + i_qpoint];
        }
        // evaluate advection operator
        (*write_face)(acc_mesh.kernel_elements());
        (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
        sw_adv.children.at("setup").stopwatch.pause();
        double dt_scaled = std::abs(basis.node(i_node) - .5)*2;
        for (int i = 0; i < 2; ++i) {
          sw_adv.children.at("BCs").stopwatch.start();
          auto& bc_cons {acc_mesh.boundary_connections()};
          #pragma omp parallel for
          for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
            int bc_sn = bc_cons[i_con].bound_cond_serial_n();
            acc_mesh.boundary_condition(bc_sn).flow_bc->apply_advection(bc_cons[i_con]);
          }
          sw_adv.children.at("BCs").stopwatch.pause();
          sw_adv.children.at("BCs").work_units_completed += acc_mesh.elements().size();
          Kernel_args ka {
            params.n_dim,
            params.row_size,
            basis,
            sw_adv.children.at("cartesian"),
            sw_adv.children.at("deformed" ),
            stopwatch.children.at("prolong/restrict"),
            acc_mesh.cartesian().kernel_connections(),
            acc_mesh.deformed ().kernel_connections(),
            acc_mesh.cartesian().kernel_elements(),
            acc_mesh.deformed ().kernel_elements(),
            acc_mesh.refined_faces(),
            dt_scaled,
            i,
            false,
            false,
          };
          compute_advection(ka);
        }
        sw_adv.children.at("cartesian").work_units_completed += acc_mesh.cartesian().elements().size();
        sw_adv.children.at("deformed" ).work_units_completed += acc_mesh.deformed ().elements().size();
        // update advection state and residual
        sw_adv.children.at("update").stopwatch.start();
        #pragma omp parallel for reduction(+:diff, n_avg)
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* state = elements[i_elem].stage(0);
          double* adv = elements[i_elem].advection_state();
          double* tss = elements[i_elem].time_step_scale();
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
            // compute update
            double pseudotime_scale = tss[i_qpoint]*2/advect_length;
            double old = adv[i_node*nq + i_qpoint]; // record for measuring residual
            adv[i_node*nq + i_qpoint] = (state[nd*nq + i_qpoint] + pseudotime_scale*1.)/(1. + pseudotime_scale);
            // add to residual
            double d = adv[i_node*nq + i_qpoint] - old;
            diff += d*d/dt_scaled/dt_scaled;
            ++n_avg;
          }
        }
        sw_adv.children.at("update").stopwatch.pause();
      }
      sw_adv.children.at("setup").work_units_completed += elements.size();
      sw_adv.children.at("update").work_units_completed += elements.size();
    }
  }
  stopwatch.children.at("set art visc").children.at("advection").stopwatch.pause();
  stopwatch.children.at("set art visc").children.at("advection").work_units_completed += elements.size();
  // update iteration status for printing to screen
  status.adv_res = std::sqrt(diff/n_avg);
  _namespace->assign("av_advection_residual", std::sqrt(diff/n_avg));
  // compute projection onto Legendre polynomial
  Eigen::VectorXd weights = basis.node_weights();
  Eigen::VectorXd orth = basis.orthogonal(av_rs - 1);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* forcing = elements[i_elem].art_visc_forcing();
    double* adv = elements[i_elem].advection_state();
    double* rk_ref = elements[i_elem].stage(1);
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double proj = 0;
      for (int i_proj = 0; i_proj < rs; ++i_proj) {
        proj += adv[i_proj*nq + i_qpoint]*weights(i_proj)*orth(i_proj);
      }
      forcing[i_qpoint] = proj*proj*2*rk_ref[(nd + 1)*nq + i_qpoint]/rk_ref[nd*nq + i_qpoint];
    }
  } // Cauchy-Kovalevskaya-style derivative estimate complete!

  // begin root-smear-square operation
  int n_real = params.n_forcing - 1; // number of real time steps (as apposed to pseudotime steps)
  double diff_time = _namespace->lookup<double>("av_diff_ratio").value()*advect_length*advect_length/n_real; // compute size of real time step (as opposed to pseudotime)
  diffuse_art_visc(n_real, diff_time);

  // clean up
  double mult = _namespace->lookup<double>("av_visc_mult").value()*advect_length;
  double us_max = advect_length*_namespace->lookup<double>("av_unscaled_max").value()*std::sqrt(2*_namespace->lookup<double>("freestream" + std::to_string(nd + 1)).value()/_namespace->lookup<double>("freestream" + std::to_string(nd)).value());
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].stage(0);
    double* rk_ref = elements[i_elem].stage(1);
    double* av = elements[i_elem].art_visc_coef();
    double* forcing = elements[i_elem].art_visc_forcing();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double f = mult*forcing[n_real*nq + i_qpoint];
      av[i_qpoint] = us_max*f/(us_max + f);
      // put the flow state back how we found it
      for (int i_var = 0; i_var < params.n_var; ++i_var) {
        int i = i_var*nq + i_qpoint;
        state[i] = rk_ref[i];
      }
    }
  }
  // update the face state
  (*write_face)(acc_mesh.kernel_elements());
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
  stopwatch.children.at("set art visc").stopwatch.pause();
  stopwatch.children.at("set art visc").work_units_completed += elements.size();
  stopwatch.stopwatch.pause();
}

void Solver::update_art_visc_elwise(double width, bool pde_based)
{
  use_art_visc = true;
  Mass mass;
  set_uncertainty(Elem_nonsmooth(mass));
  auto& elems = acc_mesh.elements();
  double scale = width/(basis.row_size - 1)*(_namespace->lookup<double>("freestream_speed").value() + _namespace->lookup<double>("freestream_sound_speed").value());
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
      double* stage0 = elems[i_elem].stage(0);
      double* stage1 = elems[i_elem].stage(1);
      double elem_av = elems[i_elem].uncertainty;
      double* av = elems[i_elem].art_visc_coef();
      double* forcing = elems[i_elem].art_visc_forcing();
      for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) stage1[i_dof] = stage0[i_dof];
      for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
        forcing[i_qpoint] = elem_av;
        forcing[params.n_qpoint() + i_qpoint] = av[i_qpoint];
      }
    }
    diffuse_art_visc(1, _namespace->lookup<double>("elementwise_art_visc_diff_ratio").value()*width*width);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      double* stage0 = elems[i_elem].stage(0);
      double* stage1 = elems[i_elem].stage(1);
      double* av = elems[i_elem].art_visc_coef();
      double* forcing = elems[i_elem].art_visc_forcing();
      for (int i_dof = 0; i_dof < params.n_dof(); ++i_dof) stage0[i_dof] = stage1[i_dof];
      for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) av[i_qpoint] = forcing[params.n_qpoint() + i_qpoint];
    }
    (*write_face)(acc_mesh.kernel_elements());
    (*kernel_factory<Prolong_refined>(params.n_dim, params.row_size, basis))(acc_mesh.refined_faces());
  } else {
    share_vertex_data([](Element& elem, int){return elem.uncertainty;},
                      [](Element& elem, int i_vert)->double&{return elem.vertex_elwise_av(i_vert);},
                      Vertex::vector_max);
    Mat<dyn, dyn> interp = Gauss_lobatto(2).interpolate(basis.nodes());
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      Eigen::Map<Mat<>> qpoint_av(elems[i_elem].art_visc_coef(), params.n_qpoint());
      Eigen::Map<Mat<>> vert_av(&elems[i_elem].vertex_elwise_av(0), params.n_vertices());
      qpoint_av = math::hypercube_matvec(interp, vert_av);
    }
  }
}

void Solver::set_art_visc_row_size(int row_size)
{
  HEXED_ASSERT(row_size >= 2, "`row_size` must be >= 2");
  HEXED_ASSERT(row_size <= basis.row_size, "`row_size` must be <= discretization row size");
  av_rs = row_size;
}

void Solver::set_fix_admissibility(bool value)
{
  fix_admis = value;
}

void Solver::set_uncertainty(const Element_func& func)
{
  auto& elems = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].uncertainty = func(elems[i_elem], basis, status.flow_time)[0];
  }
}

void Solver::set_uncert_surface_rep(int bc_sn)
{
  const int nv = params.n_var;
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int nfq = nq/params.row_size;
  // record flow state in reference stage
  // so that state storage can be used for normal vectors
  auto& elems = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].uncertainty = 0;
    elems[i_elem].record = 0;
    double* state = elems[i_elem].stage(0);
    double* ref = elems[i_elem].stage(1);
    for (int i_var = 0; i_var < nv; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        ref[i_var*nq + i_qpoint] = state[i_var*nq + i_qpoint];
      }
    }
  }
  // identify boundary elements
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    if (con.bound_cond_serial_n() == bc_sn) {
      // the `2*n_dim` identifies that this is a boundary element and the rest identifies which face is on the boundary
      // assumes each element has at most face on *this* boundary (might have other faces on other boundaries)
      con.element().record = 2*nd + 2*con.i_dim() + con.inside_face_sign();
    }
  }
  // first write the (non-unit) normal vectors to the state storage
  // and then extrapolate those to the face storage
  auto& def_elems = acc_mesh.deformed().elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < def_elems.size(); ++i_elem) {
    auto& elem = def_elems[i_elem];
    if (elem.record/(2*nd) == 1) { // if this element is on the boundary...
      // pick out the reference level normal vector which is pointing out of the flow domain
      double* state = elem.stage(0);
      double* nrml = elem.reference_level_normals() + (elem.record - 2*nd)/2*nd*nq;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
          state[i_dim*nq + i_qpoint] = nrml[i_dim*nq + i_qpoint]*math::sign(elem.record%2);
        }
      }
    }
  }
  // extrapolate to faces
  (*write_face)(acc_mesh.kernel_elements());
  // for each face, compute unit surface normal from the reference level normal
  Mat<dyn, dyn> bound = basis.boundary();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < def_elems.size(); ++i_elem) {
    auto& elem = def_elems[i_elem];
    if (elem.record/(2*nd) == 1) {
      int i_dim = (elem.record - 2*nd)/2;
      int positive = elem.record%2;
      for (int i_face = 0; i_face < 2*nd; ++i_face) {
        if (i_face/2 != i_dim) {
          // extrapolate the reference level normal to the wall surface
          // because only at the wall surface is the reference level normal the same as the wall normal
          // and then write that back to the whole face
          Eigen::Map<Mat<dyn, dyn>> face(elem.faces[i_face], nfq, nd);
          int i_dim_extrap = (nd == 3) ? i_dim > 3 - i_dim - i_face/2 : 0;
          for (int j_dim = 0; j_dim < nd; ++j_dim) {
            Mat<dyn, dyn> b = Mat<>::Ones(params.row_size)*bound(positive, all);
            Mat<> extrap = math::dimension_matvec(b, face(all, j_dim), i_dim_extrap);
            face(all, j_dim) = extrap;
          }
          // normalize to get *unit* normals
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
            face(i_qpoint, all).normalize();
          }
        }
      }
    }
  }
  (*kernel_factory<Prolong_refined>(nd, params.row_size, basis))(acc_mesh.refined_faces());
  // compute difference between neighboring elements
  auto& def_cons = acc_mesh.deformed().face_connections();
  #pragma omp parallel for
  for (int i_con = 0; i_con < def_cons.size(); ++i_con) {
    auto& con = def_cons[i_con];
    auto permute = kernel_factory<Face_permutation>(nd, params.row_size, con.direction(), con.state() + (nd + 2)*nfq);
    permute->match_faces();
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        double* f [2];
        for (int i_side : {0, 1}) f[i_side] = con.state() + (i_side*(nd + 2) + i_dim)*nfq + i_qpoint;
        *f[0] = *f[1] = *f[0] - *f[1];
      }
    }
    permute->restore();
  }
  // set difference to zero for faces that are on other boundaries
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    auto& con = bc_cons[i_con];
    double* state = con.inside_face();
    for (int i_dof = 0; i_dof < nv*nfq; ++i_dof) state[i_dof] = 0;
  }
  (*kernel_factory<Restrict_refined>(nd, params.row_size, basis, false))(acc_mesh.refined_faces());
  // total up differences for each boundary element and write to uncertainty
  // also restore the flow state to what it was at the start of this function
  Mat<> weights = math::pow_outer(basis.node_weights(), nd - 1);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    if (elem.record/(2*nd) == 1) {
      for (int i_face = 0; i_face < 2*nd; ++i_face) {
        if (i_face/2 != (elem.record - 2*nd)/2) {
          Eigen::Map<Mat<dyn, dyn>> face(elem.faces[i_face], nfq, nd);
          Mat<> norm = face.rowwise().norm();
          elem.uncertainty += std::sqrt(norm.dot(weights.asDiagonal()*norm));
        }
      }
      elem.uncertainty /= 2*std::max(nd - 1, 1);
    }
    double* state = elem.stage(0);
    double* ref = elem.stage(1);
    for (int i_var = 0; i_var < nv; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        state[i_var*nq + i_qpoint] = ref[i_var*nq + i_qpoint];
      }
    }
  }
  (*write_face)(acc_mesh.kernel_elements());
}

void Solver::synch_extruded_uncert()
{
  auto cons = acc_mesh.extruded_connections();
  bool changed = true;
  while(changed) {
    changed = false;
    #pragma omp parallel for reduction(||:changed)
    for (int i_con = 0; i_con < cons.size(); ++i_con) {
      double uncert [2];
      for (int i_side = 0; i_side < 2; ++i_side) {
        #pragma omp atomic read
        uncert[i_side] = cons[i_con].element(i_side).uncertainty;
      }
      if (uncert[0] > uncert[1] + 1e-14) {
        #pragma omp atomic write
        cons[i_con].element(1).uncertainty = uncert[0];
        changed = true;
      }
    }
  }
}

void Solver::update()
{
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  auto& elems = acc_mesh.elements();

  for (int i_flow = 0; i_flow < _namespace->lookup<int>("flow_iters").value(); ++i_flow)
  {
    // compute time step
    double safety = _namespace->lookup<double>("max_safety").value();
    double n_cheby = _namespace->lookup<double>("n_cheby_flow").value();
    double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1);
    // run chebyshev iterations
    for (int i_cheby = 0; i_cheby < n_cheby; ++i_cheby)
    {
      double nominal_dt = std::min(max_dt(safety/max_cheby, safety), _namespace->lookup<double>("max_time_step").value());
      double dt = nominal_dt*math::chebyshev_step(n_cheby, i_cheby);
      // record reference state for residual calculation
      bool fixed = false;
      // compute inviscid update
      for (int i = 0; i < 2; ++i) {
        Kernel_args ka {
          params.n_dim,
          params.row_size,
          basis,
          stopwatch.children.at("cartesian"),
          stopwatch.children.at("deformed" ),
          stopwatch.children.at("prolong/restrict"),
          acc_mesh.cartesian().kernel_connections(),
          acc_mesh.deformed ().kernel_connections(),
          acc_mesh.cartesian().kernel_elements(),
          acc_mesh.deformed ().kernel_elements(),
          acc_mesh.refined_faces(),
          dt,
          i,
          false,
          bool(_namespace->lookup<int>("use_filter").value()),
        };
        apply_state_bcs();
        if (use_ldg()) compute_viscous(dt, i, false);
        else compute_euler(ka);
        // note that function call must come first to ensure it is evaluated despite short-circuiting
        fixed = fix_admissibility(_namespace->lookup<double>("fix_admis_max_safety").value()) || fixed;
      }
      if (fixed) break;

      // update status for reporting
      _namespace->assign<double>("time_step", dt);
      _namespace->assign<double>("flow_time", _namespace->lookup<double>("flow_time").value() + dt);
      status.time_step = dt;
      status.flow_time += dt;
    }
  }

  _namespace->assign("iteration", _namespace->lookup<int>("iteration").value() + 1);
  _namespace->assign("wall_time", status.wall_time());
  ++status.iteration;
  stopwatch.work_units_completed += elems.size();
  stopwatch.children.at("cartesian").work_units_completed += acc_mesh.cartesian().elements().size();
  stopwatch.children.at("deformed" ).work_units_completed += acc_mesh.deformed ().elements().size();
  stopwatch.stopwatch.pause();
}

void Solver::update_implicit()
{
  HEXED_ASSERT(_implicit, "`update_implicit` called on a Solver that was not constructed in implicit mode");
  Linearized lin(*this);
  iterative::gmres(lin, 27, 1);
  lin.add(-Linearized::storage_start, 1., -Linearized::storage_start + 3, 1., 0.);
  (*write_face)(acc_mesh.kernel_elements());
  (*kernel_factory<Prolong_refined>(params.n_dim, params.row_size, basis))(acc_mesh.refined_faces());
  fix_admissibility(.7);
}

void Solver::compute_residual()
{
  apply_state_bcs();
  Kernel_args ka {
    params.n_dim,
    params.row_size,
    basis,
    stopwatch.children.at("cartesian"),
    stopwatch.children.at("deformed" ),
    stopwatch.children.at("prolong/restrict"),
    acc_mesh.cartesian().kernel_connections(),
    acc_mesh.deformed ().kernel_connections(),
    acc_mesh.cartesian().kernel_elements(),
    acc_mesh.deformed ().kernel_elements(),
    acc_mesh.refined_faces(),
    1.,
    0,
    true,
    bool(_namespace->lookup<int>("use_filter").value()),
  };
  if (use_ldg()) compute_viscous(1., 0, true);
  else compute_euler(ka);
}

Iteration_status Solver::iteration_status()
{
  Iteration_status stat = status;
  return stat;
}

bool Solver::is_admissible()
{
  auto& sw = stopwatch.children.at("check admis.");
  sw.stopwatch.start();
  auto& elems = acc_mesh.elements();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  bool admiss = 1;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].record = 0;
  }
  #pragma omp parallel for reduction (&&:admiss)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    bool elem_admis = true;
    elem_admis = elem_admis && hexed::thermo::admissible(elem.stage(0), nd, nq);
    for (int i_face = 0; i_face < params.n_dim*2; ++i_face) {
      elem_admis = elem_admis && hexed::thermo::admissible(elem.faces[i_face], nd, nq/rs);
    }
    if (!elem_admis) elem.record = 1;
    admiss = admiss && elem_admis;
  }
  auto& ref_faces = acc_mesh.refined_faces();
  bool refined_admiss = 1;
  #pragma omp parallel for reduction (&&:refined_admiss)
  for (int i_face = 0; i_face < ref_faces.size(); ++i_face) {
    auto& ref = ref_faces[i_face];
    int n_fine = params.n_vertices()/2;
    for (int i_dim = 0; i_dim < nd - 1; ++i_dim) n_fine /= 1 + ref.stretch[i_dim];
    for (int i_fine = 0; i_fine < n_fine; ++i_fine) {
      refined_admiss = refined_admiss && hexed::thermo::admissible(ref.fine[i_fine], nd, nq/rs);
    }
  }
  sw.work_units_completed += acc_mesh.elements().size();
  sw.stopwatch.pause();
  return admiss && refined_admiss;
}

bool Solver::fix_admissibility(double stability_ratio)
{
  if (!fix_admis) return false;
  auto& sw_fix = stopwatch.children.at("fix admis.");
  sw_fix.stopwatch.start();
  std::string wd = _namespace->lookup<std::string>("working_dir").value();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  const int nv = params.n_vertices();
  int iter;
  int n_iters = std::numeric_limits<int>::max();
  for (iter = 0; iter < n_iters; ++iter) {
    HEXED_ASSERT(iter < 1e5, format_str(200, "failed to fix thermodynamic admissability in %i iterations", iter));
    #if HEXED_USE_XDMF
    if (iter == 100) {
      State_variables sv;
      Record rec;
      std::vector<const Qpoint_func*> to_vis {&sv, &rec};
      visualize_field("xdmf", wd + "severe_indamis" + std::to_string(status.iteration), Qf_concat(to_vis));
    }
    #endif
    if (is_admissible()) {
      if (iter) n_iters = std::min(n_iters, 2*iter);
      else {
        ++iter;
        break;
      }
    } else {
      n_iters = std::numeric_limits<int>::max();
    }
    if (iter == 0) {
      _printer->print(format_str(200, "Thermodynamically inadmissible state detected (solver iteration %i). Attempting to fix...\n",
                                 _namespace->lookup<int>("iteration").value()));
    }
    auto bounds = bounds_field(State_variables(), 2*rs);
    _printer->print(format_str(200, "    iteration %i: mass in [%e, %e]; energy in [%e, %e]\n", iter, bounds[nd][0], bounds[nd][1], bounds[nd + 1][0], bounds[nd + 1][1]));
    auto& elems = acc_mesh.elements();
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_vert = 0; i_vert < nv; ++i_vert) {
        elem.vertex_fix_admis_coef(i_vert) = elem.record;
      }
    }
    share_vertex_data(&Element::vertex_fix_admis_coef, Vertex::vector_max);
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
    share_vertex_data(&Element::vertex_fix_admis_coef, Vertex::vector_max);
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
      Eigen::Map<Mat<>>(elem.fix_admis_coef(), nq) = math::hypercube_matvec(interp, vert_fac);
    }
    #if HEXED_USE_XDMF
    if (status.iteration >= last_fix_vis_iter + 1000 && iter == 0) {
      last_fix_vis_iter = status.iteration;
      State_variables sv;
      Record rec;
      Fix_admis_coef fac;
      std::vector<const Qpoint_func*> to_vis {&sv, &rec, &fac};
      visualize_field("xdmf", wd + "inadmis" + std::to_string(status.iteration), Qf_concat(to_vis));
    }
    #endif
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        std::swap(elem.fix_admis_coef()[i_qpoint], elem.art_visc_coef()[i_qpoint]);
      }
    }
    double dt = stability_ratio;
    (*kernel_factory<Spatial<pde::Fix_therm_admis, false>::Max_dt>(nd, rs, basis, true, false, dt, dt))(acc_mesh.cartesian().kernel_elements(), stopwatch.children.at("fix admis.").children.at("cartesian"), "compute time step");
    (*kernel_factory<Spatial<pde::Fix_therm_admis,  true>::Max_dt>(nd, rs, basis, true, false, dt, dt))(acc_mesh.deformed ().kernel_elements(), stopwatch.children.at("fix admis.").children.at("deformed" ), "compute time step");
    dt = 1.;
    double linear = dt;
    double quadratic = dt*dt/8/0.9;
    std::array<double, 2> step;
    step[1] = (linear + std::sqrt(linear*linear - 4*quadratic))/2.;
    step[0] = quadratic/step[1];
    for (double s : step) {
      auto& bc_cons {acc_mesh.boundary_connections()};
      #pragma omp parallel for
      for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
        double* in_f = bc_cons[i_con].inside_face();
        double* gh_f = bc_cons[i_con].ghost_face();
        for (int i_dof = 0; i_dof < nq*(nd + 2)/rs; ++i_dof) gh_f[i_dof] = in_f[i_dof];
      }
      compute_fta(s, 0);
    }
    double safety = _namespace->lookup<double>("max_safety").value();
    double n_cheby = _namespace->lookup<double>("n_cheby_flow").value();
    double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1);
    max_dt(safety/max_cheby, safety);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      auto& elem = elems[i_elem];
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        std::swap(elem.fix_admis_coef()[i_qpoint], elem.art_visc_coef()[i_qpoint]);
      }
    }
  }
  --iter;
  if (iter) _printer->print("done\n");
  status.fix_admis_iters += iter;
  _namespace->assign("fix_iters", _namespace->lookup<int>("fix_iters").value() + iter);
  sw_fix.work_units_completed += acc_mesh.elements().size()*iter;
  sw_fix.stopwatch.pause();
  return iter;
}

void Solver::reset_counters()
{
  status.fix_admis_iters = 0;
  _namespace->assign("fix_iters", 0);
}

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func& func)
{
  return func(acc_mesh.element(ref_level, is_deformed, serial_n), basis, i_qpoint, status.flow_time);
}

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, const Element_func& func)
{
  return func(acc_mesh.element(ref_level, is_deformed, serial_n), basis, status.flow_time);
}

std::vector<double> Solver::integral_field(const Qpoint_func& integrand)
{
  // compute `n_dim`-dimensional quadrature weights from 1D weights
  Eigen::VectorXd weights = math::pow_outer(basis.node_weights(), params.n_dim);
  // now compute the integral with the above quadrature weights
  std::vector<double> integral (integrand.n_var(params.n_dim), 0.);
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    Element& element {elements[i_elem]};
    double volume = math::pow(element.nominal_size(), params.n_dim);
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
      auto qpoint_integrand {integrand(element, basis, i_qpoint, status.flow_time)};
      for (unsigned i_var = 0; i_var < qpoint_integrand.size(); ++i_var) {
        integral[i_var] += weights[i_qpoint]*volume*qpoint_integrand[i_var]*element.jacobian_determinant(i_qpoint);
      }
    }
  }
  return integral;
}

std::vector<double> Solver::integral_surface(const Boundary_func& integrand, int bc_sn)
{
  // setup
  const int nd = params.n_dim;
  const int n_int = integrand.n_var(nd);
  const int nq = params.n_qpoint();
  const int nfq = nq/basis.row_size;
  Eigen::MatrixXd boundary = basis.boundary();
  Eigen::VectorXd weights = math::pow_outer(basis.node_weights(), params.n_dim - 1);
  // write the state to the faces so that the BCs can access it
  (*write_face)(acc_mesh.kernel_elements());
  // compute the integral
  std::vector<double> integral (n_int, 0.);
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con)
  {
    auto& con {bc_cons[i_con]};
    auto& elem = con.element();
    double area = math::pow(elem.nominal_size(), nd - 1);
    if (con.bound_cond_serial_n() == bc_sn)
    {
      double* nrml = con.normal();
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        double nrml_mag = 0;
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          nrml_mag += math::pow(nrml[i_dim*nfq + i_qpoint], 2);
        }
        nrml_mag = std::sqrt(nrml_mag);
        auto qpoint_integrand = integrand(con, i_qpoint, status.flow_time);
        for (int i_var = 0; i_var < n_int; ++i_var) {
          integral[i_var] += qpoint_integrand[i_var]*weights(i_qpoint)*area*nrml_mag;
        }
      }
    }
  }
  return integral;
};

std::vector<std::array<double, 2>> Solver::bounds_field(const Qpoint_func& func, int n_sample)
{
  const int n_var = func.n_var(params.n_dim);
  std::vector<std::array<double, 2>> bounds(n_var);
  for (int i_var = 0; i_var < n_var; ++i_var) {
    bounds[i_var] = {std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
  }
  const int n_block = math::pow(n_sample, params.n_dim);
  auto& elems = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem)
  {
    Element& elem {elems[i_elem]};
    Eigen::VectorXd vars = Vis_data(elem, func, basis, status.flow_time).interior(n_sample);
    for (int i_var = 0; i_var < n_var; ++i_var) {
      auto var = vars(Eigen::seqN(i_var*n_block, n_block));
      bounds[i_var][0] = std::min(var.minCoeff(), bounds[i_var][0]);
      bounds[i_var][1] = std::max(var.maxCoeff(), bounds[i_var][1]);
    }
  }
  return bounds;
}

std::unique_ptr<Visualizer> Solver::_visualizer(std::string format, std::string name, const Output_data& output_variables, int n_dim_topo)
{
  std::unique_ptr<Visualizer> visualizer;
  if (format == "xdmf") {
    #if HEXED_USE_XDMF
    visualizer.reset(new Xdmf_wrapper(params.n_dim, n_dim_topo, name, output_variables, status.flow_time));
    #else
    HEXED_ASSERT(false, "`format = xdmf` requires `USE_XDMF ON`");
    #endif
  } else if (format == "tecplot") {
    #if HEXED_USE_TECPLOT
    const int n_vis = output_variables.n_var(params.n_dim); // number of variables to visualize
    std::vector<std::string> var_names;
    for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back(output_variables.variable_name(params.n_dim, i_vis));
    visualizer.reset(new Tecplot_file(name, params.n_dim, n_dim_topo, var_names, status.flow_time));
    #else
    HEXED_ASSERT(false, "`format = tecplot` requires `USE_TECPLOT ON`");
    #endif
  } else HEXED_ASSERT(false, format_str(1000, "visualization format `%s` not recognized", format.c_str()));
  return visualizer;
}

void Solver::visualize_field(std::string format, std::string name, const Qpoint_func& output_variables, int n_sample, bool wireframe)
{
  auto visualizer = _visualizer(format, name, output_variables, wireframe ? 1 : params.n_dim);
  HEXED_ASSERT(params.n_dim > wireframe, "can only visualize field wireframes in > 1D");
  Position_func pos_func;
  int nv = output_variables.n_var(params.n_dim);
  auto& elems = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Vis_data pos_dat(elems[i_elem], pos_func, basis, status.flow_time);
    Vis_data out_dat(elems[i_elem], output_variables, basis, status.flow_time);
    if (wireframe) {
      Mat<> pos = pos_dat.edges(n_sample);
      Mat<> out = out_dat.edges(n_sample);
      for (int i_edge = 0; i_edge < math::pow(2, params.n_dim - 1)*params.n_dim; ++i_edge) {
        visualizer->write_block(n_sample, pos.data() + i_edge*params.n_dim*n_sample, out.data() + i_edge*nv*n_sample);
      }
    } else {
      Mat<> pos = pos_dat.interior(n_sample);
      Mat<> out = out_dat.interior(n_sample);
      visualizer->write_block(n_sample, pos.data(), out.data());
    }
  }
}

void Solver::visualize_surface(std::string format, std::string name, int bc_sn, const Boundary_func& func, int n_sample, bool wireframe)
{
  HEXED_ASSERT(params.n_dim > 1, "cannot visualize surfaces in 1D");
  HEXED_ASSERT(params.n_dim > 1 + wireframe, "can only visualize surface wireframes in 3D");
  auto visualizer = _visualizer(format, name, func, wireframe ? 1 : params.n_dim - 1);
  // convenience definitions
  const int nfq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = func.n_var(nd);
  const int n_block {math::pow(n_sample, nd - 1)};
  const int n_edge {math::pow(n_sample, nd - 2)};
  // setup
  Mat<dyn, dyn> interp = basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.));
  Mat<dyn, dyn> boundary = basis.boundary();
  // iterate through boundary connections and visualize a zone for each
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con)
  {
    auto& con {bc_cons[i_con]};
    if (con.bound_cond_serial_n() == bc_sn)
    {
      auto& elem = con.element();
      // fetch the position
      const int i_face = con.direction().i_face(0);
      Mat<dyn, dyn> qpoint_pos (nfq, nd);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        auto pos = elem.face_position(basis, i_face, i_qpoint);
        for (int i_dim = 0; i_dim < nd; ++i_dim) qpoint_pos(i_qpoint, i_dim) = pos[i_dim];
      }
      // fetch the output variables
      Mat<dyn, dyn> qpoint_vars (nfq, nv);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        auto vars = func(con, i_qpoint, status.flow_time);
        for (int i_var = 0; i_var < nv; ++i_var) {
          qpoint_vars(i_qpoint, i_var) = vars[i_var];
        }
      }
      if (wireframe) {
        for (int i_dim = 0; i_dim < params.n_dim - 1; ++i_dim) {
          for (int i_sign = 0; i_sign < 2; ++i_sign) {
            Mat<dyn, dyn> interp_pos (n_edge, nd);
            for (int j_dim = 0; j_dim < nd; ++j_dim) {
              // extrapolate from faces to edges
              Mat<> uniform = math::dimension_matvec(boundary(i_sign, all), qpoint_pos.col(j_dim), i_dim);
              // interpolate from edge qpoints to uniformly-spaced sample points
              interp_pos.col(j_dim) = math::hypercube_matvec(interp, uniform);
            }
            Mat<dyn, dyn> interp_vars (n_edge, nv);
            for (int i_var = 0; i_var < nv; ++i_var) {
              Mat<> uniform = math::dimension_matvec(boundary(i_sign, all), qpoint_vars.col(i_var), i_dim);
              interp_vars.col(i_var) = math::hypercube_matvec(interp, uniform);
            }
            visualizer->write_block(n_sample, interp_pos.data(), interp_vars.data());
          }
        }
      } else {
        // interpolate from quadrature points to sample points
        Mat<dyn, dyn> interp_pos (n_block, nd);
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          interp_pos.col(i_dim) = math::hypercube_matvec(interp, qpoint_pos.col(i_dim));
        }
        Mat<dyn, dyn> interp_vars (n_block, nv);
        for (int i_var = 0; i_var < nv; ++i_var) {
          interp_vars.col(i_var) = math::hypercube_matvec(interp, qpoint_vars.col(i_var));
        }
        // visualize
        visualizer->write_block(n_sample, interp_pos.data(), interp_vars.data());
      }
    }
  }
  visualizer.reset();
}

void Solver::vis_cart_surf(std::string format, std::string name, int bc_sn, const Boundary_func& func)
{
  Mesh::Reset_vertices reset(acc_mesh);
  visualize_surface(format, name, bc_sn, func, 2);
}

void Solver::vis_lts_constraints(std::string format, std::string name, int n_sample)
{
  auto& elems = acc_mesh.elements();
  int nf = params.n_dof();
  int nq = params.n_qpoint();
  double n_cheby = _namespace->lookup<double>("n_cheby_flow").value();
  double max_cheby = math::chebyshev_step(n_cheby, n_cheby - 1);
  // write local time steps for convection and diffusion to the mass and energy of the reference state.
  // Reference state is used for storage because `Element::time_step_scale` only has space for one scalar
  for (int i_term = 0; i_term < 2; ++i_term) {
    double safeties [] {1., huge};
    max_dt(safeties[i_term]/max_cheby, safeties[!i_term]);
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
      Eigen::Map<Mat<>>(elems[i_elem].stage(1) + (params.n_dim + i_term)*nq, nq) = Eigen::Map<Mat<>>(elems[i_elem].time_step_scale(), nq);
    }
  }
  // swap current state and reference state
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Eigen::Map<Mat<>> stage0(elems[i_elem].stage(0), nf);
    Eigen::Map<Mat<>> stage1(elems[i_elem].stage(1), nf);
    Mat<> ref_state = stage1;
    stage1 = stage0;
    stage0 = ref_state;
  }
  // visualize. Note that visualizing straight from the reference state would require implementing another `Qpoint_func` which would be ugly
  Interpreter inter(std::vector<std::string>{});
  Struct_expr expr("lts_convective = density; lts_diffusive = energy; lts_ratio = lts_diffusive/lts_convective;");
  visualize_field(format, name, Qpoint_expr(expr, inter), n_sample);
  // restore the current state from the reference state
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    Eigen::Map<Mat<>>(elems[i_elem].stage(0), nf) = Eigen::Map<Mat<>>(elems[i_elem].stage(1), nf);
  }
}

}
