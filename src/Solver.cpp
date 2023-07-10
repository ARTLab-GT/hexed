#include <fstream>
#include <iostream>

#include <config.hpp>
#include <Solver.hpp>
#include <Tecplot_file.hpp>
#include <Vis_data.hpp>
#include <otter_vis.hpp>
#include <thermo.hpp>

// kernels
#include <Restrict_refined.hpp>
#include <Prolong_refined.hpp>
#include <Face_permutation.hpp>

#include <pde.hpp>
#include <Spatial.hpp>

namespace hexed
{

void Solver::share_vertex_data(Element::vertex_value_access access_func, Vertex::reduction reduce)
{
  auto& elements = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].push_shareable_value(access_func);
  }
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].fetch_shareable_value(access_func, reduce);
  }
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  #pragma omp parallel for
  for (int i_match = 0; i_match < matchers.size(); ++i_match) matchers[i_match].match(access_func);
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

double Solver::max_dt()
{
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  if (use_ldg()) {
    return std::min((*kernel_factory<Spatial<Element         , pde::Navier_stokes<true >::Pde>::Max_dt>(nd, rs, basis, is_local_time, visc, therm_cond))(acc_mesh.cartesian().elements(), sw_car, "compute time step"),
                    (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<true >::Pde>::Max_dt>(nd, rs, basis, is_local_time, visc, therm_cond))(acc_mesh.deformed ().elements(), sw_def, "compute time step"));
  } else {
    return std::min((*kernel_factory<Spatial<Element         , pde::Navier_stokes<false>::Pde>::Max_dt>(nd, rs, basis, is_local_time                  ))(acc_mesh.cartesian().elements(), sw_car, "compute time step"),
                    (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<false>::Pde>::Max_dt>(nd, rs, basis, is_local_time                  ))(acc_mesh.deformed ().elements(), sw_def, "compute time step"));
  }
}

Solver::Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping,
               Transport_model viscosity_model, Transport_model thermal_conductivity_model) :
  params{3, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size},
  basis{row_size},
  stopwatch{"(element*iteration)"},
  use_art_visc{false},
  fix_admis{false},
  av_rs{basis.row_size},
  write_face(kernel_factory<Spatial<Element, pde::Navier_stokes<false>::Pde>::Write_face>(params.n_dim, params.row_size, basis)), // note: for now, false and true are equivalent for `Write_face`
  is_local_time{local_time_stepping},
  visc{viscosity_model},
  therm_cond{thermal_conductivity_model}
{
  // setup categories for performance reporting
  stopwatch.children.emplace("initialize reference", stopwatch.work_unit_name);
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
  stopwatch.children.emplace("residual computation", "element*evaluation");
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

void Solver::relax_vertices(double factor)
{
  // relax the vertices
  auto verts {acc_mesh.vertices()};
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].calc_relax(factor);
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].apply_relax();
}

void Solver::snap_vertices()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).mesh_bc->snap_vertices(bc_cons[i_con]);
  }
  // if any immobile vertices have strayed from their nominal position (probably by `eat`ing)
  // snap them back where they belong
  auto& elems = acc_mesh.cartesian().elements();
  int nd = params.n_dim;
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      Lock::Acquire a(elem.vertex(i_vert).lock);
      auto& pos = elem.vertex(i_vert).pos;
      double nom_sz = elem.nominal_size();
      auto nom_pos = elem.nominal_position();
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        pos[i_dim] = nom_sz*(nom_pos[i_dim] + (i_vert/math::pow(2, nd - 1 - i_dim))%2);
      }
      pos(Eigen::seqN(0, nd)) += elem.origin;
    }
  }
  // vertex relaxation/snapping will cause hanging vertices to drift away from hanging vertex faces they are supposed to be coincident with
  // so now we put them back where they belong
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  #pragma omp parallel for
  for (int i_match = 0; i_match < matchers.size(); ++i_match) {
    matchers[i_match].match(&Element::vertex_position<0>);
    matchers[i_match].match(&Element::vertex_position<1>);
    matchers[i_match].match(&Element::vertex_position<2>);
  }
}

void Solver::snap_faces()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    Lock::Acquire acq(bc_cons[i_con].element().lock);
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).mesh_bc->snap_node_adj(bc_cons[i_con], basis);
  }
}

void Solver::calc_jacobian()
{
  acc_mesh.valid().assert_valid();
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
  if (func.n_var(params.n_dim) != params.n_var) {
    throw std::runtime_error("initializer has wrong number of output variables");
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
  (*write_face)(elements);
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

void Solver::set_art_visc_smoothness(double advect_length)
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
      double scale = std::sqrt(2*rk_ref[(nd + 1)*nq + i_qpoint]*rk_ref[nd*nq + i_qpoint]);
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
  (*write_face)(elements);
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
  (*kernel_factory<Spatial<Element         , pde::Advection>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), sw_adv.children.at("cartesian"), "compute time step");
  (*kernel_factory<Spatial<Deformed_element, pde::Advection>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), sw_adv.children.at("deformed" ), "compute time step");
  double dt_adv = av_advect_stab_rat;

  // begin estimation of high-order derivative in the style of the Cauchy-Kovalevskaya theorem using a linear advection equation.
  double diff = 0; // for residual computation
  int n_avg = 0;
  // perform pseudotime iteration
  for (int iter = 0; iter < av_advect_iters; ++iter)
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
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
            state[nd*nq + i_qpoint] = state[(nd + 1)*nq + i_qpoint] = adv[i_node*nq + i_qpoint];
          }
        }
        // evaluate advection operator
        (*write_face)(elements);
        (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
        sw_adv.children.at("setup").stopwatch.pause();
        double dt_scaled = dt_adv*std::abs(basis.node(i_node) - .5)*2;
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
          compute_advection(dt_scaled, i);
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
            double pseudotime_scale = + dt_adv*tss[i_qpoint]*2/advect_length;
            double old = adv[i_node*nq + i_qpoint]; // record for measuring residual
            adv[i_node*nq + i_qpoint] += state[nd*nq + i_qpoint] - state[(nd + 1)*nq + i_qpoint] + pseudotime_scale*1.;
            adv[i_node*nq + i_qpoint] /= 1. + pseudotime_scale;
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
  // compute projection onto Legendre polynomial
  Eigen::VectorXd weights = basis.node_weights();
  Eigen::VectorXd orth = basis.orthogonal(av_rs - 1);
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* forcing = elements[i_elem].art_visc_forcing();
    double* adv = elements[i_elem].advection_state();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      double proj = 0;
      for (int i_proj = 0; i_proj < rs; ++i_proj) {
        proj += adv[i_proj*nq + i_qpoint]*weights(i_proj)*orth(i_proj);
      }
      forcing[i_qpoint] = proj*proj;
    }
  } // Cauchy-Kovalevskaya-style derivative estimate complete!

  // begin root-smear-square operation
  int n_real = 3; // number of real time steps (as apposed to pseudotime steps)
  // evaluate CFL condition
  (*kernel_factory<Spatial<Element         , pde::Smooth_art_visc>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("cartesian"), "compute time step");
  (*kernel_factory<Spatial<Deformed_element, pde::Smooth_art_visc>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("deformed" ), "compute time step");
  double dt_diff = av_diff_stab_rat;
  double diff_time = av_diff_ratio/n_real; // compute size of real time step (as opposed to pseudotime)
  // initialize residual to zero (will compute RMS over all real time steps)
  status.diff_res = 0;
  for (int real_step = 0; real_step < n_real; ++real_step)
  {
    // copy value from forcing function storage to scalar variable of stage 0
    #pragma omp parallel for
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      double* state = elements[i_elem].stage(0);
      double* forcing = elements[i_elem].art_visc_forcing();
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        state[i_qpoint] = forcing[(real_step + 1)*nq + i_qpoint];
      }
    }
    (*write_face)(elements);
    (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
    // set up multistage scheme
    double linear = dt_diff;
    double quadratic = basis.cancellation_diffusive()/basis.max_cfl_diffusive()*dt_diff*dt_diff;
    std::array<double, 2> step;
    step[1] = (linear + std::sqrt(linear*linear - 4*quadratic))/2.;
    step[0] = quadratic/step[1];
    diff = 0;
    n_avg = 0;
    // perform pseudotime iteration
    for (int i_iter = 0; i_iter < av_diff_iters; ++i_iter)
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
      for (double s : step)
      {
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
            double pseudotime_scale = s/diff_time*tss[i_qpoint]/advect_length;
            state[i_qpoint] += forcing[real_step*nq + i_qpoint]*pseudotime_scale;
            state[i_qpoint] /= 1 + pseudotime_scale;
          }
        }
        // update face state
        (*write_face)(elements);
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
        diff += d*d/dt_diff/dt_diff;
        ++n_avg;
        // update forcing
        forcing[(real_step + 1)*nq + i_qpoint] = state[i_qpoint];
      }
    }
    status.diff_res += diff/n_avg;
  }
  status.diff_res = std::sqrt(status.diff_res/n_real); // finish computing RMS residual
  // clean up
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    double* state = elements[i_elem].stage(0);
    double* rk_ref = elements[i_elem].stage(1);
    double* av = elements[i_elem].art_visc_coef();
    double* forcing = elements[i_elem].art_visc_forcing();
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      // set artificial viscosity to square root of diffused scalar state times scaling factor
      double scale_sq = 2*rk_ref[(nd + 1)*nq + i_qpoint]/rk_ref[nd*nq + i_qpoint];
      double f = std::max(0., forcing[n_real*nq + i_qpoint]);
      f = av_unscaled_max*f/(av_unscaled_max + f);
      av[i_qpoint] = av_visc_mult*advect_length*std::sqrt(f*scale_sq); // root-smear-square complete!
      // put the flow state back how we found it
      for (int i_var = 0; i_var < params.n_var; ++i_var) {
        int i = i_var*nq + i_qpoint;
        state[i] = rk_ref[i];
      }
    }
  }
  // update the face state
  (*write_face)(elements);
  (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
  stopwatch.children.at("set art visc").stopwatch.pause();
  stopwatch.children.at("set art visc").work_units_completed += elements.size();
  stopwatch.stopwatch.pause();
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

void Solver::set_resolution_badness(const Element_func& func)
{
  auto& elems = acc_mesh.elements();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].resolution_badness = func(elems[i_elem], basis, status.flow_time)[0];
  }
}

void Solver::set_res_bad_surface_rep(int bc_sn)
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
    elems[i_elem].resolution_badness = 0;
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
  (*write_face)(elems);
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
  // total up differences for each boundary element and write to resolution badness
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
          elem.resolution_badness += std::sqrt(norm.dot(weights.asDiagonal()*norm));
        }
      }
      elem.resolution_badness /= 2*std::max(nd - 1, 1);
    }
    double* state = elem.stage(0);
    double* ref = elem.stage(1);
    for (int i_var = 0; i_var < nv; ++i_var) {
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        state[i_var*nq + i_qpoint] = ref[i_var*nq + i_qpoint];
      }
    }
  }
  (*write_face)(elems);
}

void Solver::synch_extruded_res_bad()
{
  auto cons = acc_mesh.extruded_connections();
  bool changed = true;
  while(changed) {
    changed = false;
    #pragma omp parallel for reduction(||:changed)
    for (int i_con = 0; i_con < cons.size(); ++i_con) {
      double bad [2];
      for (int i_side = 0; i_side < 2; ++i_side) {
        #pragma omp atomic read
        bad[i_side] = cons[i_con].element(i_side).resolution_badness;
      }
      if (bad[0] > bad[1] + 1e-14) {
        #pragma omp atomic write
        cons[i_con].element(1).resolution_badness = bad[0];
        changed = true;
      }
    }
  }
}

void Solver::update(double dt, bool cfl_driven)
{
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  auto& elems = acc_mesh.elements();

  // compute time step
  double max = max_dt();
  if (cfl_driven) dt *= max;
  else if (dt > max) std::cerr << "warning: time step exceeds CFL limit\n";

  // record reference state for Runge-Kutta scheme
  const int n_dof = params.n_dof();
  auto& irk = stopwatch.children.at("initialize reference");
  irk.stopwatch.start();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].stage(0);
    for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] = state[i_dof];
  }
  irk.stopwatch.pause();
  irk.work_units_completed += elems.size();

  // compute inviscid update
  for (int i = 0; i < 2; ++i) {
    apply_state_bcs();
    if (use_ldg()) compute_viscous(dt, i);
    else compute_inviscid(dt, i);
    fix_admissibility(fix_admis_stab_rat);
  }

  // update status for reporting
  status.time_step = dt;
  status.flow_time += dt;
  ++status.iteration;
  stopwatch.stopwatch.pause();
  stopwatch.work_units_completed += elems.size();
  stopwatch.children.at("cartesian").work_units_completed += acc_mesh.cartesian().elements().size();
  stopwatch.children.at("deformed" ).work_units_completed += acc_mesh.deformed ().elements().size();
}

Iteration_status Solver::iteration_status()
{
  Iteration_status stat = status;
  stopwatch.children.at("residual computation").stopwatch.start();
  Physical_update update;
  auto res = integral_field(Pow(update, 2));
  for (double& r : res) r /= stat.time_step*stat.time_step;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    stat.mmtm_res += res[i_dim];
  }
  stat.mmtm_res = std::sqrt(stat.mmtm_res);
  stat.mass_res = std::sqrt(res[params.n_dim]);
  stat.ener_res = std::sqrt(res[params.n_dim + 1]);
  stopwatch.children.at("residual computation").stopwatch.pause();
  stopwatch.children.at("residual computation").work_units_completed += acc_mesh.elements().size();
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
  #pragma omp parallel for reduction (&&:admiss)
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    auto& elem = elems[i_elem];
    admiss = admiss && hexed::thermo::admissible(elem.stage(0), nd, nq);
    for (int i_face = 0; i_face < params.n_dim*2; ++i_face) {
      admiss = admiss && hexed::thermo::admissible(elem.faces[i_face], nd, nq/rs);
    }
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

void Solver::fix_admissibility(double stability_ratio)
{
  if (!fix_admis) return;
  auto& sw_fix = stopwatch.children.at("fix admis.");
  sw_fix.stopwatch.start();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint();
  const int rs = params.row_size;
  int iter;
  for (iter = 0;; ++iter) {
    if (iter > 99999) {
      #if HEXED_USE_OTTER
      otter::plot plt;
      visualize_field_otter(plt, Pressure(), 1, {0, 0}, Pressure(), {0, 0}, otter::const_colormap(Eigen::Vector4d{1., 0., 0., .1}), otter::plasma, false, false);
      visualize_field_otter(plt, Pressure(), 0);
      visualize_edges_otter(plt);
      plt.show();
      #endif
      char buffer [200];
      snprintf(buffer, 200, "failed to fix thermodynamic admissability in %i iterations", iter);
      throw std::runtime_error(buffer);
    }
    if (is_admissible()) break;
    else {
      if (status.iteration >= last_fix_vis_iter + 1000) {
        last_fix_vis_iter = status.iteration;
      }
      (*kernel_factory<Spatial<Element         , pde::Fix_therm_admis>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), stopwatch.children.at("fix admis.").children.at("cartesian"), "compute time step");
      (*kernel_factory<Spatial<Deformed_element, pde::Fix_therm_admis>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), stopwatch.children.at("fix admis.").children.at("deformed" ), "compute time step");
      double dt = stability_ratio;
      double linear = dt;
      double quadratic = basis.cancellation_diffusive()/basis.max_cfl_diffusive()*dt*dt;
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
      max_dt();
    }
  }
  status.fix_admis_iters += iter;
  sw_fix.work_units_completed += acc_mesh.elements().size()*iter;
  sw_fix.stopwatch.pause();
}

void Solver::reset_counters()
{
  status.fix_admis_iters = 0;
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
  (*write_face)(acc_mesh.elements());
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

#if HEXED_USE_TECPLOT
void Solver::visualize_field_tecplot(const Qpoint_func& output_variables, std::string name, int n_sample, bool include_edges, bool include_qpoints, bool include_interior)
{
  const int n_dim = params.n_dim;
  const int n_vis = output_variables.n_var(n_dim); // number of variables to visualize
  const int n_corners {math::pow(2, n_dim - 1)};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  std::vector<std::string> var_names;
  for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back(output_variables.variable_name(n_dim, i_vis));
  Tecplot_file file {name, n_dim, var_names, status.flow_time};

  auto& elems = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem)
  {
    Element& elem {elems[i_elem]};
    Vis_data vis_pos(elem, hexed::Position_func(), basis, status.flow_time);
    Vis_data vis_out(elem, output_variables, basis, status.flow_time);
    // note: each visualization stage is enclosed in `{}` to ensure that only one `Tecplot_file::Zone` is alive at a time
    // visualize edges
    if (n_dim > 1 && include_edges) { // 1D elements don't really have edges
      Tecplot_file::Line_segments edges {file, n_dim*n_corners, n_sample, "edges"};
      auto edge_pos = vis_pos.edges(n_sample);
      auto edge_state = vis_out.edges(n_sample);
      for (int i_edge = 0; i_edge < n_corners*n_dim; ++i_edge) {
        edges.write(edge_pos.data() + i_edge*n_sample*n_dim, edge_state.data() + i_edge*n_sample*n_vis);
      }
    }

    if (include_qpoints) { // visualize quadrature points
      Tecplot_file::Structured_block qpoints {file, basis.row_size, "element_qpoints"};
      qpoints.write(vis_pos.qpoints().data(), vis_out.qpoints().data());
    }

    if (include_interior) { // visualize interior (that is, quadrature point data interpolated to a fine mesh of sample points)
      Tecplot_file::Structured_block interior {file, n_sample, "element_interior"};
      auto interp_pos = vis_pos.interior(n_sample);
      auto interp_out = vis_out.interior(n_sample);
      interior.write(interp_pos.data(), interp_out.data());
    }
  }
}

void Solver::visualize_field_tecplot(std::string name, int n_sample, bool edges, bool qpoints, bool interior)
{
  State_variables sv;
  Art_visc_coef avc;
  Advection_state as(av_rs);
  Art_visc_forcing avf;
  std::vector<const Qpoint_func*> funcs {&sv};
  if (use_art_visc) {
    funcs.push_back(&avc);
    funcs.push_back(&as);
    funcs.push_back(&avf);
  }
  visualize_field_tecplot(Qf_concat(funcs), name, n_sample, edges, qpoints, interior);
}

void Solver::visualize_surface_tecplot(int bc_sn, const Boundary_func& func, std::string name, int n_sample)
{
  if (params.n_dim == 1) throw std::runtime_error("cannot visualize surfaces in 1D");
  // convenience definitions
  const int nfq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = func.n_var(nd);
  const int n_block {math::pow(n_sample, nd - 1)};
  // setup
  std::vector<std::string> var_names;
  for (int i_var = 0; i_var < nv; ++i_var) var_names.push_back(func.variable_name(nd, i_var));
  Tecplot_file file {name, nd, var_names, status.flow_time};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
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
      // interpolate from quadrature points to sample points
      Mat<dyn, dyn> interp_pos (n_block, nd);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        interp_pos.col(i_dim) = math::hypercube_matvec(interp, qpoint_pos.col(i_dim));
      }
      // fetch the state
      Mat<dyn, dyn> qpoint_vars (nfq, nv);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        auto vars = func(con, i_qpoint, status.flow_time);
        for (int i_var = 0; i_var < nv; ++i_var) {
          qpoint_vars(i_qpoint, i_var) = vars[i_var];
        }
      }
      // interpolate to sample points
      Mat<dyn, dyn> interp_vars (n_block, nv);
      for (int i_var = 0; i_var < nv; ++i_var) {
        interp_vars.col(i_var) = math::hypercube_matvec(interp, qpoint_vars.col(i_var));
      }
      // visualize zone
      Tecplot_file::Structured_block zone {file, n_sample, "face_interior", nd - 1};
      zone.write(interp_pos.data(), interp_vars.data());
    }
  }
  status.flow_time += 1; // FIXME
}

void Solver::visualize_surface_tecplot(int bc_sn, std::string name, int n_sample)
{
  State_variables sv;
  Outward_normal on;
  Viscous_stress vs;
  Heat_flux hf;
  std::vector<const Boundary_func*> funcs {&sv, &on};
  if (visc.is_viscous) funcs.push_back(&vs);
  if (therm_cond.is_viscous) funcs.push_back(&hf);
  visualize_surface_tecplot(bc_sn, Bf_concat(funcs), name, n_sample);
}

void Solver::vis_cart_surf_tecplot(int bc_sn, std::string name)
{
  auto verts = acc_mesh.vertices();
  Mat<dyn, dyn> pos(3, verts.size());
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    pos(all, i_vert) = verts[i_vert].pos;
  }
  acc_mesh.reset_vertices();
  visualize_surface_tecplot(bc_sn, hexed::Resolution_badness(), name, 2);
  #pragma omp parallel for
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) {
    verts[i_vert].pos = pos(all, i_vert);
  }
}
#endif

#if HEXED_USE_OTTER
void Solver::visualize_edges_otter(otter::plot& plt, Eigen::Matrix<double, 1, Eigen::Dynamic> color, int n_sample)
{
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    otter_vis::add_edges(plt, elements[i_elem], basis, color, n_sample);
  }
}

void Solver::visualize_surface_otter(otter::plot& plt, int bc_sn, const otter::colormap& cmap, const Qpoint_func& color_by, std::array<double, 2> bounds, bool transparent, int n_sample, double tol)
{
  if (color_by.n_var(params.n_dim) != 1) throw std::runtime_error("`color_by` must be scalar");
  // substitute bounds if necessary
  if (std::isnan(bounds[0]) || std::isnan(bounds[1])) {
    bounds = bounds_field(color_by)[0];
    bounds[1] += tol;
  }
  // iterate through boundary connections and visualize an `otter::surface` for each
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con)
  {
    auto& con {bc_cons[i_con]};
    if (con.bound_cond_serial_n() == bc_sn)
    {
      // fetch face data
      Vis_data vis_pos(con.element(), Position_func(), basis, status.flow_time);
      Vis_data vis_var(con.element(),   color_by, basis, status.flow_time);
      // interpolate to boundary face
      auto dir = con.direction();
      Eigen::MatrixXd face_pos = vis_pos.face(dir.i_dim[0], dir.face_sign[0], n_sample);
      Eigen::VectorXd face_var = vis_var.face(dir.i_dim[0], dir.face_sign[0], n_sample);
      // transform so that bounds[0] is the bottom of the color map and bounds[1] is the top
      face_var = (face_var - Eigen::VectorXd::Constant(face_var.size(), bounds[0]))/(bounds[1] - bounds[0]);
      // add to plot
      if (params.n_dim == 2) {
        face_pos.resize(n_sample, params.n_dim);
        otter::curve curve(face_pos, cmap(face_var));
        plt.add(curve);
      } else if (params.n_dim == 3) {
        face_pos.resize(n_sample*n_sample, params.n_dim);
        Eigen::MatrixXd color = Eigen::MatrixXd::Constant(face_var.size(), 4, transparent ? .2 : 1.);
        Eigen::MatrixXd mapped = cmap(face_var);
        color(Eigen::all, Eigen::seqN(0, mapped.cols())) = mapped;
        otter::surface surf(n_sample, face_pos, color);
        if (transparent) surf.transparency = otter::surface::transparent_add;
        plt.add(surf);
      }
    }
  }
}

void Solver::visualize_field_otter(otter::plot& plt,
                                   const Qpoint_func& contour,
                                   int n_contour,
                                   std::array<double, 2> contour_bounds,
                                   const Qpoint_func& color_by,
                                   std::array<double, 2> color_bounds,
                                   const otter::colormap& cmap_contour,
                                   const otter::colormap& cmap_field,
                                   bool transparent, bool show_color,
                                   int n_div, double tol)
{
  if (contour.n_var(params.n_dim) != 1) throw std::runtime_error("`contour` must be scalar");
  if (color_by.n_var(params.n_dim) != 1) throw std::runtime_error("`color_by` must be scalar");
  // substitute bounds if necessary
  auto bf = bounds_field(contour)[0];
  if (std::isnan(contour_bounds[0]) || std::isnan(contour_bounds[1])) {
    contour_bounds = bf;
    contour_bounds[1] += tol;
  }
  if (std::isnan(color_bounds[0]) || std::isnan(color_bounds[1])) {
    color_bounds = bounds_field(color_by)[0];
    color_bounds[1] += tol;
  }
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    auto& elem = elements[i_elem];
    for (int i_contour = 0; i_contour < n_contour; ++i_contour) {
      double contour_val = (i_contour + 1)/(n_contour + 1.);
      contour_val = contour_val*(contour_bounds[1] - contour_bounds[0]) + contour_bounds[0];
      contour_val = (contour_val - bf[0])/(bf[1] - bf[0]);
      // add contour line/surface
      otter_vis::add_contour(plt, elem, basis, Scaled(contour, bf), contour_val, n_div,
                             color_by, color_bounds, cmap_contour, transparent, status.flow_time, tol);
      // for 2d, color the flow field as well
    }
    if ((params.n_dim == 2) && show_color) {
      const int n_sample = 2*n_div + 1;
      Vis_data vis_pos(elem, Position_func(), basis, status.flow_time);
      Vis_data vis_var(elem,         contour, basis, status.flow_time);
      Eigen::MatrixXd pos_3d(n_sample*n_sample, 3);
      Eigen::MatrixXd pos = vis_pos.interior(n_sample);
      Eigen::MatrixXd var = vis_var.interior(n_sample);
      pos_3d(Eigen::all, Eigen::seqN(0, 2)) = pos;
      // set z componento to slightly negative so that lines show up on top
      pos_3d(Eigen::all, 2).array() = -1e-4*elem.nominal_size();
      // adjust color-by variable so that bound interval maps to [0, 1]
      var = (var - Eigen::VectorXd::Constant(var.size(), contour_bounds[0]))/(contour_bounds[1] - contour_bounds[0]);
      // throw it up there
      plt.add(otter::surface(n_sample, pos_3d, cmap_field(var)));
    }
  }
}
#endif

}
