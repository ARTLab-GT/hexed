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
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].push_shareable_value(access_func);
  }
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].fetch_shareable_value(access_func, reduce);
  }
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  for (int i_match = 0; i_match < matchers.size(); ++i_match) matchers[i_match].match(access_func);
}

void Solver::apply_state_bcs()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).flow_bc->apply_state(bc_cons[i_con]);
  }
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

void Solver::apply_avc_diff_flux_bcs()
{
  int nd = params.n_dim;
  int rs = params.row_size;
  int nq = params.n_qpoint();
  auto& bc_cons {acc_mesh.boundary_connections()};
  #pragma omp parallel for
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* in_f = bc_cons[i_con].inside_face() + 2*(nd + 2)*nq/rs;
    double* gh_f = bc_cons[i_con].ghost_face()  + 2*(nd + 2)*nq/rs;
    for (int i_qpoint = 0; i_qpoint < nq/rs; ++i_qpoint) gh_f[i_qpoint] = -in_f[i_qpoint];
  }
}

bool Solver::use_ldg()
{
  return visc || use_art_visc;
}

double Solver::max_dt()
{
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  if (use_ldg()) {
    return std::min((*kernel_factory<Spatial<Element         , pde::Navier_stokes<true >::Pde>::Max_dt>(nd, rs, basis, is_local_time, visc))(acc_mesh.cartesian().elements(), sw_car, "compute time step"),
                    (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<true >::Pde>::Max_dt>(nd, rs, basis, is_local_time, visc))(acc_mesh.deformed ().elements(), sw_def, "compute time step"));
  } else {
    return std::min((*kernel_factory<Spatial<Element         , pde::Navier_stokes<false>::Pde>::Max_dt>(nd, rs, basis, is_local_time      ))(acc_mesh.cartesian().elements(), sw_car, "compute time step"),
                    (*kernel_factory<Spatial<Deformed_element, pde::Navier_stokes<false>::Pde>::Max_dt>(nd, rs, basis, is_local_time      ))(acc_mesh.deformed ().elements(), sw_def, "compute time step"));
  }
}

Solver::Solver(int n_dim, int row_size, double root_mesh_size, bool local_time_stepping, bool viscous) :
  params{3, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size},
  basis{row_size},
  stopwatch{"(element*iteration)"},
  use_art_visc{false},
  fix_admis{false},
  av_rs{basis.row_size},
  write_face(kernel_factory<Spatial<Element, pde::Navier_stokes<false>::Pde>::Write_face>(params.n_dim, params.row_size, basis)), // note: for now, false and true are equivalent for `Write_face`
  visc{viscous},
  is_local_time{local_time_stepping}
{
  // setup categories for performance reporting
  stopwatch.children.emplace("initialize reference", stopwatch.work_unit_name);
  std::string unit = "(element*(time integration stage))";
  stopwatch.children.emplace("prolong/restrict", unit);
  stopwatch.children.emplace("fix admis.", "(element*(fix admis. iter))");
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
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].calc_relax(factor);
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].apply_relax();
}

void Solver::snap_vertices()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    int bc_sn = bc_cons[i_con].bound_cond_serial_n();
    acc_mesh.boundary_condition(bc_sn).mesh_bc->snap_vertices(bc_cons[i_con]);
  }
  // vertex relaxation/snapping will cause hanging vertices to drift away from hanging vertex faces they are supposed to be coincident with
  // so now we put them back where they belong
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  for (int i_match = 0; i_match < matchers.size(); ++i_match) {
    matchers[i_match].match(&Element::vertex_position<0>);
    matchers[i_match].match(&Element::vertex_position<1>);
    matchers[i_match].match(&Element::vertex_position<2>);
  }
}

void Solver::snap_faces()
{
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
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
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].set_jacobian(basis);
  }

  /*
   * compute surface normals for deformed connections
   */
  auto& def_cons {acc_mesh.deformed().face_connections()};
  for (int i_con = 0; i_con < def_cons.size(); ++i_con) {
    double* nrml = def_cons[i_con].normal();
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) nrml[i_data] = 0.;
  }
  // for deformed refined faces, set normal to coarse face normal (for Cartesian, setting normal is not necessary)
  auto& ref_cons = acc_mesh.deformed().refined_connections();
  (*kernel_factory<Prolong_refined>(n_dim, rs, basis, true))(acc_mesh.refined_faces());
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
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* in_f = bc_cons[i_con].inside_face();
    double* gh_f = bc_cons[i_con].ghost_face();
    for (int i_data = 0; i_data < n_dim*nfq; ++i_data) gh_f[i_data] = in_f[i_data];
  }
  // compute the shared face normal
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
}

void Solver::set_art_visc_off()
{
  use_art_visc = false;
}

void Solver::set_art_visc_constant(double value)
{
  use_art_visc = true;
  auto& elements = acc_mesh.elements();
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
  double heat_rat = 1.4;
  const int nq = params.n_qpoint();
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& elements = acc_mesh.elements();
  auto& bc_cons {acc_mesh.boundary_connections()};

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
      double scale = 2*heat_rat*rk_ref[(nd + 1)*nq + i_qpoint]*rk_ref[nd*nq + i_qpoint];
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        double mmtm = rk_ref[i_dim*nq + i_qpoint];
        scale += (1. - heat_rat)*mmtm*mmtm;
      }
      scale = std::sqrt(scale);
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
  (*kernel_factory<Spatial<Element         , pde::Advection>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), sw_adv.children.at("cartesian"), "max char speed");
  (*kernel_factory<Spatial<Deformed_element, pde::Advection>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), sw_adv.children.at("deformed" ), "max char speed");
  double dt_adv = av_advect_stab_rat;
  // ensure that time step is small enough that the time derivative term will be stable
  {
    double max_rat = 0;
    #pragma omp parallel for reduction(max:max_rat)
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      double* tss = elements[i_elem].time_step_scale();
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        max_rat = std::max(max_rat, dt_adv/advect_length*tss[i_qpoint]);
      }
    }
    dt_adv /= std::max(max_rat/av_advect_stab_rat, 1.);
  }

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
        State_variables sv;
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
            double d = state[nd*nq + i_qpoint] - state[(nd + 1)*nq + i_qpoint]
                       + dt_adv*tss[i_qpoint]*(1. - adv[i_node*nq + i_qpoint])*2/advect_length;
            adv[i_node*nq + i_qpoint] += d; // add to advection state
            // add to residual
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
  (*kernel_factory<Spatial<Element         , pde::Smooth_art_visc>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("cartesian"), "max char speed"),
  (*kernel_factory<Spatial<Deformed_element, pde::Smooth_art_visc>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), stopwatch.children.at("set art visc").children.at("diffusion").children.at("deformed" ), "max char speed");
  double dt_diff = av_diff_stab_rat;
  double diff_time = av_diff_ratio*advect_length/n_real; // compute size of real time step (as opposed to pseudotime)
  // adjust for time derivative term if necessary
  {
    double max_rat = 0;
    #pragma omp parallel for reduction(max:max_rat)
    for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
      double* tss = elements[i_elem].time_step_scale();
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        max_rat = std::max(max_rat, dt_diff/diff_time*tss[i_qpoint]*tss[i_qpoint]);
      }
    }
    dt_diff /= std::max(max_rat/av_diff_stab_rat, 1.);
  }
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
        // compute (real) time derivative term
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* state = elements[i_elem].stage(0);
          double* forcing = elements[i_elem].art_visc_forcing();
          double* tss = elements[i_elem].time_step_scale();
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
            state[nq + i_qpoint] = s*(forcing[real_step*nq + i_qpoint] - state[i_qpoint])/diff_time*tss[i_qpoint]*tss[i_qpoint];
          }
        }
        // apply value boundary condition (or lack thereof, in this case)
        #pragma omp parallel for
        for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
          double* in_f = bc_cons[i_con].inside_face();
          double* gh_f = bc_cons[i_con].ghost_face();
          for (int i_dof = 0; i_dof < nq/rs; ++i_dof) {
            gh_f[i_dof] = in_f[i_dof];
          }
        }
        compute_avc_diff(s, 0);
        // update state
        #pragma omp parallel for
        for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
          double* state = elements[i_elem].stage(0);
          for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
            state[i_qpoint] += state[nq + i_qpoint];
          }
        }
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
      // set artificial viscosity to square root of diffused scalar state times stagnation enthalpy (with user-defined multiplier)
      double mass = rk_ref[nd*nq + i_qpoint];
      double scale_sq = 2*heat_rat*rk_ref[(nd + 1)*nq + i_qpoint]/mass;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        double veloc = rk_ref[i_dim*nq + i_qpoint]/mass;
        scale_sq += (1. - heat_rat)*veloc*veloc;
      }
      av[i_qpoint] = av_visc_mult*advect_length*std::sqrt(std::max(0., forcing[n_real*nq + i_qpoint]*scale_sq)); // root-smear-square complete!
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
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    elems[i_elem].resolution_badness = func(elems[i_elem], basis, status.flow_time)[0];
  }
}

void Solver::synch_extruded_res_bad()
{
  auto cons = acc_mesh.extruded_connections();
  bool changed = true;
  while(changed) {
    changed = false;
    for (int i_con = 0; i_con < cons.size(); ++i_con) {
      double* bad [2];
      for (int i_side = 0; i_side < 2; ++i_side) {
        bad[i_side] = &cons[i_con].element(i_side).resolution_badness;
      }
      if (*bad[0] > *bad[1] + 1e-14) {
        *bad[1] = *bad[0];
        changed = true;
      }
    }
  }
}

void Solver::update(double stability_ratio)
{
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  auto& elems = acc_mesh.elements();

  // compute time step
  double dt = stability_ratio*max_dt();

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
    fix_admissibility(fix_admis_stab_rat*stability_ratio);
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

void Solver::fix_admissibility(double stability_ratio)
{
  if (!fix_admis) return;
  auto& sw_fix = stopwatch.children.at("fix admis.");
  sw_fix.stopwatch.start();
  auto& elems = acc_mesh.elements();
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
    if (admiss && refined_admiss) break;
    else {
      (*kernel_factory<Spatial<Element         , pde::Fix_therm_admis>::Max_dt>(nd, rs, basis, true))(acc_mesh.cartesian().elements(), stopwatch.children.at("fix admis.").children.at("cartesian"), "max char speed");
      (*kernel_factory<Spatial<Deformed_element, pde::Fix_therm_admis>::Max_dt>(nd, rs, basis, true))(acc_mesh.deformed ().elements(), stopwatch.children.at("fix admis.").children.at("deformed" ), "max char speed");
      double dt = stability_ratio;
      double linear = dt;
      double quadratic = basis.cancellation_diffusive()/basis.max_cfl_diffusive()*dt*dt;
      std::array<double, 2> step;
      step[1] = (linear + std::sqrt(linear*linear - 4*quadratic))/2.;
      step[0] = quadratic/step[1];
      for (double s : step) {
        apply_state_bcs();
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
  Eigen::VectorXd weights = custom_math::pow_outer(basis.node_weights(), params.n_dim);
  // now compute the integral with the above quadrature weights
  std::vector<double> integral (integrand.n_var(params.n_dim), 0.);
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    Element& element {elements[i_elem]};
    double volume = custom_math::pow(element.nominal_size(), params.n_dim);
    for (int i_qpoint = 0; i_qpoint < params.n_qpoint(); ++i_qpoint) {
      auto qpoint_integrand {integrand(element, basis, i_qpoint, status.flow_time)};
      for (unsigned i_var = 0; i_var < qpoint_integrand.size(); ++i_var) {
        integral[i_var] += weights[i_qpoint]*volume*qpoint_integrand[i_var]*element.jacobian_determinant(i_qpoint);
      }
    }
  }
  return integral;
}

std::vector<double> Solver::integral_surface(const Surface_func& integrand, int bc_sn)
{
  // setup
  const int nd = params.n_dim;
  const int nv = params.n_var;
  const int n_int = integrand.n_var(nd);
  const int nq = params.n_qpoint();
  const int nfq = nq/basis.row_size;
  Eigen::MatrixXd boundary = basis.boundary();
  Eigen::VectorXd weights = custom_math::pow_outer(basis.node_weights(), params.n_dim - 1);
  // write the state to the faces so that the BCs can access it
  (*write_face)(acc_mesh.elements());
  // compute the integral
  std::vector<double> integral (n_int, 0.);
  auto& bc_cons {acc_mesh.boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con)
  {
    auto& con {bc_cons[i_con]};
    auto& elem = con.element();
    double area = custom_math::pow(elem.nominal_size(), nd - 1);
    if (con.bound_cond_serial_n() == bc_sn)
    {
      double* state = con.state();
      double* nrml = con.normal();
      double* pos = con.surface_position();
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint)
      {
        // fetch state
        std::vector<double> qpoint_state;
        for (int i_var = 0; i_var < nv; ++i_var) {
          qpoint_state.push_back(state[i_var*nfq + i_qpoint]);
        }
        // fetch position and normal vector
        std::vector<double> qpoint_pos;
        std::vector<double> normal;
        double nrml_mag = 0;
        int nrml_sign = 2*con.direction().flip_normal(0) - 1;
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          qpoint_pos.push_back(pos[i_dim*nfq + i_qpoint]);
          double comp = nrml_sign*nrml[i_dim*nfq + i_qpoint];
          normal.push_back(comp);
          nrml_mag += comp*comp;
        }
        nrml_mag = std::sqrt(nrml_mag);
        for (int i_dim = 0; i_dim < nd; ++i_dim) normal[i_dim] /= nrml_mag;
        auto qpoint_integrand = integrand(qpoint_pos, status.flow_time, qpoint_state, normal);
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
  const int n_block = custom_math::pow(n_sample, params.n_dim);
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
void Solver::visualize_field_tecplot(const Qpoint_func& output_variables, std::string name, int n_sample)
{
  const int n_dim = params.n_dim;
  const int n_vis = output_variables.n_var(n_dim); // number of variables to visualize
  const int n_corners {custom_math::pow(2, n_dim - 1)};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  std::vector<std::string> var_names;
  for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back(output_variables.variable_name(i_vis));
  Tecplot_file file {name, n_dim, var_names, status.flow_time};

  for (int i_elem = 0; i_elem < acc_mesh.elements().size(); ++i_elem)
  {
    Element& elem {acc_mesh.elements()[i_elem]};
    Vis_data vis_pos(elem, hexed::Position_func(), basis, status.flow_time);
    Vis_data vis_out(elem, output_variables, basis, status.flow_time);
    // note: each visualization stage is enclosed in `{}` to ensure that only one `Tecplot_file::Zone` is alive at a time
    // visualize edges
    if (n_dim > 1) // 1D elements don't really have edges
    {
      Tecplot_file::Line_segments edges {file, n_dim*n_corners, n_sample, "edges"};
      auto edge_pos = vis_pos.edges(n_sample);
      auto edge_state = vis_out.edges(n_sample);
      for (int i_edge = 0; i_edge < n_corners*n_dim; ++i_edge) {
        edges.write(edge_pos.data() + i_edge*n_sample*n_dim, edge_state.data() + i_edge*n_sample*n_vis);
      }
    }

    { // visualize quadrature points
      Tecplot_file::Structured_block qpoints {file, basis.row_size, "element_qpoints"};
      qpoints.write(vis_pos.qpoints().data(), vis_out.qpoints().data());
    }

    { // visualize interior (that is, quadrature point data interpolated to a fine mesh of sample points)
      Tecplot_file::Structured_block interior {file, n_sample, "element_interior"};
      auto interp_pos = vis_pos.interior(n_sample);
      auto interp_out = vis_out.interior(n_sample);
      interior.write(interp_pos.data(), interp_out.data());
    }
  }
}

void Solver::visualize_surface_tecplot(int bc_sn, std::string name, int n_sample)
{
  if (params.n_dim == 1) throw std::runtime_error("cannot visualize surfaces in 1D");
  // convenience definitions
  const int nfq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = params.n_var;
  const int n_block {custom_math::pow(n_sample, nd - 1)};
  // setup
  std::vector<std::string> var_names;
  for (int i_var = 0; i_var < nv; ++i_var) var_names.push_back(State_variables().variable_name(i_var));
  Tecplot_file file {name, nd, var_names, status.flow_time};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  // write the state to the faces so that the BCs can access it
  (*write_face)(acc_mesh.elements());
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
      Eigen::MatrixXd qpoint_pos (nfq, nd);
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        auto pos = elem.face_position(basis, i_face, i_qpoint);
        for (int i_dim = 0; i_dim < nd; ++i_dim) qpoint_pos(i_qpoint, i_dim) = pos[i_dim];
      }
      // interpolate from quadrature points to sample points
      Eigen::MatrixXd interp_pos (n_block, nd);
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        interp_pos.col(i_dim) = custom_math::hypercube_matvec(interp, qpoint_pos.col(i_dim));
      }
      // fetch/interpolate the state
      Eigen::MatrixXd interp_state (n_block, nv);
      for (int i_var = 0; i_var < nv; ++i_var) {
        Eigen::Map<Eigen::VectorXd> qpoint_state (con.state() + i_var*nfq, nfq);
        interp_state.col(i_var) = custom_math::hypercube_matvec(interp, qpoint_state);
      }
      // visualize zone
      Tecplot_file::Structured_block zone {file, n_sample, "face_interior", nd - 1};
      zone.write(interp_pos.data(), interp_state.data());
    }
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
