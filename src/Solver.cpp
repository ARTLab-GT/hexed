#include <config.hpp>
#include <Solver.hpp>
#include <Mcs_cartesian.hpp>
#include <Mcs_deformed.hpp>
#include <Write_face.hpp>
#include <Prolong_refined.hpp>
#include <Neighbor_cartesian.hpp>
#include <Neighbor_deformed.hpp>
#include <Restrict_refined.hpp>
#include <Local_cartesian.hpp>
#include <Local_deformed.hpp>
#include <Face_permutation.hpp>
#include <Tecplot_file.hpp>
#include <Vis_data.hpp>
#include <otter_vis.hpp>

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

Solver::Solver(int n_dim, int row_size, double root_mesh_size) :
  params{2, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size},
  basis{row_size},
  stopwatch{"(element*iteration)"}
{
  // setup categories for performance reporting
  stopwatch.children.emplace("initialize RK", stopwatch.work_unit_name);
  std::string unit = "(element*RK stage)";
  stopwatch.children.emplace("prolong/restrict", unit);
  for (std::string type : {"cartesian", "deformed"}) {
    stopwatch.children.emplace(type, stopwatch.work_unit_name);
    auto& children = stopwatch.children.at(type).children;
    children.emplace("max char speed", stopwatch.work_unit_name);
    children.emplace("neighbor", "(connection*RK stage)");
    children.emplace("local", unit);
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
  const int n_dim = params.n_dim;
  const int n_jac = n_dim*n_dim;
  const int rs = basis.row_size;
  const int nfq = params.n_qpoint()/rs;

  // compute element jacobians
  auto& elements = acc_mesh.deformed().elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    elements[i_elem].set_jacobian(basis);
  }

  // compute Jacobian for deformed connections
  auto& def_cons {acc_mesh.deformed().face_connections()};
  for (int i_con = 0; i_con < def_cons.size(); ++i_con) {
    double* jac = def_cons[i_con].jacobian();
    for (int i_data = 0; i_data < n_jac*nfq; ++i_data) jac[i_data] = 0.;
  }
  auto bound_mat = basis.boundary();
  // compute 1 entry at a time
  for (int i_dim = 0; i_dim < n_dim; ++i_dim)
  {
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      // compute the face jacobian of each element
      for (int i_elem = 0; i_elem < elements.size(); ++i_elem)
      {
        auto& elem = elements[i_elem];
        Eigen::Map<Eigen::VectorXd> elem_jac (elem.jacobian() + (i_dim*n_dim + j_dim)*params.n_qpoint(), params.n_qpoint());
        for (int k_dim = 0; k_dim < n_dim; ++k_dim) {
          for (int i_side : {0, 1}) {
            // use the element face storage to temporarily store one entry of the extrapolated jacobian
            Eigen::Map<Eigen::VectorXd> face_jac (elem.face() + (2*k_dim + i_side)*params.n_var*nfq, nfq);
            face_jac = custom_math::dimension_matvec(bound_mat.row(i_side), elem_jac, k_dim); // extrapolate jacobian to faces
          }
        }
      }
      // prolong Jacobian onto fine faces at hanging node connections
      auto& ref_cons = acc_mesh.deformed().refined_connections();
      // scale jacobian on refined faces with stretching
      for (int i_ref = 0; i_ref < ref_cons.size(); ++i_ref) {
        auto& ref = ref_cons[i_ref];
        auto& rf = ref.refined_face;
        int face_dim = ref.direction().i_dim[ref.order_reversed()];
        // if the column of the jacobian currently being processed is one that needs to be stretched...
        if ((j_dim != face_dim) && rf.stretch[(n_dim == 3) && (2*j_dim > 3 - face_dim)]) {
          // double this component of the jacobian for the benefit of the fine faces
          double* data = rf.coarse_face();
          for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) data[i_qpoint] *= 2;
        }
      }
      (*kernel_factory<Prolong_refined>(n_dim, rs, basis))(acc_mesh.deformed().refined_faces());
      // for BCs, copy Jacobian to ghost face
      auto& def_bc_cons = acc_mesh.deformed().boundary_connections();
      for (int i_con = 0; i_con < def_bc_cons.size(); ++i_con) {
        double* in_f = def_bc_cons[i_con].inside_face();
        double* gh_f = def_bc_cons[i_con].ghost_face();
        for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) gh_f[i_qpoint] = in_f[i_qpoint];
      }
      // compute the shared face jacobian
      for (int i_con = 0; i_con < def_cons.size(); ++i_con)
      {
        auto& con = def_cons[i_con];
        double* shared_jac = con.jacobian();
        double* elem_jac [2];
        int normal_sign [2]; // record any sign change due to normal flipping
        auto dir = con.direction();
        for (int i_side : {0, 1}) {
          elem_jac[i_side] = con.face(i_side);
          normal_sign[i_side] = ((j_dim == dir.i_dim[i_side]) && (dir.flip_normal(i_side))) ? -1 : 1;
        }
        // permute face 1 so that quadrature points match up
        (*kernel_factory<Face_permutation>(n_dim, rs, dir, elem_jac[1])).match_faces();
        for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
          // take average of element face jacobians with appropriate axis permutations
          int col = j_dim;
          int tangential_sign = 1;
          // swap dimension `dir.i_dim[0]` and dimension `dir.i_dim[1]`
          if (j_dim == dir.i_dim[1]) {
            col = dir.i_dim[0];
          }
          if (j_dim == dir.i_dim[0]) {
            col = dir.i_dim[1];
            if (dir.flip_tangential()) tangential_sign = -1;
          }
          shared_jac[(i_dim*n_dim + j_dim)*nfq + i_qpoint] += 0.5*normal_sign[0]*elem_jac[0][i_qpoint];
          shared_jac[(i_dim*n_dim + col)*nfq + i_qpoint] += 0.5*normal_sign[1]*tangential_sign*elem_jac[1][i_qpoint];
        }
      }
    }
  }

  // set Jacobian of Cartesian boundary connections to identity
  auto& bc_cons {acc_mesh.cartesian().boundary_connections()};
  for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
    double* jac = bc_cons[i_con].jacobian();
    for (int i_jac = 0; i_jac < n_jac; ++i_jac) {
      int sign = ((i_jac%params.n_dim == bc_cons[i_con].i_dim()) && !bc_cons[i_con].inside_face_sign()) ? -1 : 1;
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
        jac[i_jac*nfq + i_qpoint] = sign*(i_jac/params.n_dim == i_jac%params.n_dim); // 1 if row == col else 0
      }
    }
  }
}

void Solver::set_local_tss()
{
  // set time step to be continuous at vertices
  share_vertex_data(&Element::vertex_time_step_scale, Vertex::vector_min);
  // construct a matrix for 1D linear interpolation
  Eigen::MatrixXd lin_interp {basis.row_size, 2};
  for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint) {
    double node {basis.node(i_qpoint)};
    lin_interp(i_qpoint, 0) = 1. - node;
    lin_interp(i_qpoint, 1) = node;
  }
  // interpolate tss to quadrature points
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    Element& elem = elements[i_elem];
    Eigen::VectorXd vert_tss (params.n_vertices());
    for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
      vert_tss(i_vert) = elem.vertex_time_step_scale(i_vert);
    }
    Eigen::Map<Eigen::VectorXd> qpoint_tss (elem.time_step_scale(), params.n_qpoint());
    qpoint_tss = custom_math::hypercube_matvec(lin_interp, vert_tss);
  }
}

void Solver::initialize(const Spacetime_func& func)
{
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
  (*kernel_factory<Write_face>(params.n_dim, params.row_size, basis))(elements);
}

void Solver::update(double stability_ratio)
{
  stopwatch.stopwatch.start(); // ready or not the clock is countin'
  auto& sw_car {stopwatch.children.at("cartesian")};
  auto& sw_def {stopwatch.children.at("deformed" )};
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& elems = acc_mesh.elements();
  // compute time step
  double mcs = std::max((*kernel_factory<Mcs_cartesian>(nd, rs))(acc_mesh.cartesian().elements(), sw_car, "max char speed"),
                        (*kernel_factory<Mcs_deformed >(nd, rs))(acc_mesh.deformed ().elements(), sw_def, "max char speed"));
  double dt = stability_ratio*basis.max_cfl_convective()/params.n_dim/mcs;
  // record reference state for Runge-Kutta scheme
  const int n_dof = params.n_dof();
  auto& irk = stopwatch.children.at("initialize RK");
  irk.stopwatch.start();
  #pragma omp parallel for
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].stage(0);
    for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] = state[i_dof];
  }
  irk.stopwatch.pause();
  irk.work_units_completed += elems.size();
  // perform update for each Runge-Kutta stage
  for (double rk_weight : rk_weights) {
    (*kernel_factory<Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
    auto& bc_cons {acc_mesh.boundary_connections()};
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      int bc_sn = bc_cons[i_con].bound_cond_serial_n();
      acc_mesh.boundary_condition(bc_sn).flow_bc->apply(bc_cons[i_con]);
    }
    (*kernel_factory<Neighbor_cartesian>(nd, rs))(acc_mesh.cartesian().face_connections(), sw_car, "neighbor");
    (*kernel_factory<Neighbor_deformed >(nd, rs))(acc_mesh.deformed ().face_connections(), sw_def, "neighbor");
    (*kernel_factory<Restrict_refined>(nd, rs, basis))(acc_mesh.refined_faces(), stopwatch.children.at("prolong/restrict"));
    (*kernel_factory<Local_cartesian>(nd, rs, basis, dt, rk_weight))(acc_mesh.cartesian().elements(), sw_car, "local");
    (*kernel_factory<Local_deformed >(nd, rs, basis, dt, rk_weight))(acc_mesh.deformed ().elements(), sw_def, "local");
  }
  // update status for reporting
  status.time_step = dt;
  status.flow_time += dt;
  ++status.iteration;
  stopwatch.stopwatch.pause();
  stopwatch.work_units_completed += elems.size();
  sw_car.work_units_completed += acc_mesh.cartesian().elements().size();
  sw_def.work_units_completed += acc_mesh.deformed ().elements().size();
}

std::vector<double> Solver::sample(int ref_level, bool is_deformed, int serial_n, int i_qpoint, const Qpoint_func& func)
{
  return func(acc_mesh.element(ref_level, is_deformed, serial_n), basis, i_qpoint, status.flow_time);
}

std::vector<double> Solver::integral_field(const Qpoint_func& integrand)
{
  /*
   * compute the quadrature point weights as an `n_dim`-dimensional outer product
   * I added some extra rows/columns of zeros so that `hypercube_matvec` doesn't get confused
   * about the dimensionality.
   */
  Eigen::VectorXd one {Eigen::VectorXd::Zero(params.n_vertices())};
  one(0) = 1.;
  Eigen::MatrixXd weights_1d {Eigen::MatrixXd::Zero(params.row_size, 2)};
  weights_1d.col(0) = basis.node_weights();
  Eigen::VectorXd weights = custom_math::hypercube_matvec(weights_1d, one);
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
  Eigen::VectorXd one {Eigen::VectorXd::Zero(params.n_vertices()/2)};
  one(0) = 1.;
  Eigen::MatrixXd weights_1d {Eigen::MatrixXd::Zero(params.row_size, 2)};
  weights_1d.col(0) = basis.node_weights();
  Eigen::VectorXd weights = custom_math::hypercube_matvec(weights_1d, one);
  // write the state to the faces so that the BCs can access it
  (*kernel_factory<Write_face>(nd, params.row_size, basis))(acc_mesh.elements());
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
      double* state = con.inside_face();
      double* jac = con.jacobian();
      for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint)
      {
        std::vector<double> qpoint_pos {elem.face_position(basis, con.direction().i_face(0), i_qpoint)};
        std::vector<double> qpoint_state;
        for (int i_var = 0; i_var < nv; ++i_var) {
          qpoint_state.push_back(state[i_var*nfq + i_qpoint]);
        }
        Eigen::MatrixXd qpoint_jac (nd, nd);
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          for (int j_dim = 0; j_dim < nd; ++j_dim) {
            qpoint_jac(i_dim, j_dim) = jac[(i_dim*nd + j_dim)*nfq + i_qpoint];
          }
        }
        Eigen::MatrixXd orth = custom_math::orthonormal(qpoint_jac, con.i_dim());
        std::vector<double> normal;
        for (int i_dim = 0; i_dim < nd; ++i_dim) normal.push_back(-orth(i_dim, con.i_dim()));
        auto qpoint_integrand = integrand(qpoint_pos, status.flow_time, qpoint_state, normal);
        qpoint_jac.col(con.i_dim()) = orth.col(con.i_dim());
        double jac_det = std::abs(qpoint_jac.determinant());
        for (int i_var = 0; i_var < n_int; ++i_var) {
          integral[i_var] += qpoint_integrand[i_var]*weights(i_qpoint)*area*jac_det;
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
  (*kernel_factory<Write_face>(nd, params.row_size, basis))(acc_mesh.elements());
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
        Eigen::Map<Eigen::VectorXd> qpoint_state (con.inside_face() + i_var*nfq, nfq);
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
                                   bool transparent,
                                   int n_div, double tol)
{
  if (contour.n_var(params.n_dim) != 1) throw std::runtime_error("`contour` must be scalar");
  if (color_by.n_var(params.n_dim) != 1) throw std::runtime_error("`color_by` must be scalar");
  // substitute bounds if necessary
  if (std::isnan(contour_bounds[0]) || std::isnan(contour_bounds[1])) {
    contour_bounds = bounds_field(contour)[0];
    contour_bounds[1] += tol;
  }
  if (std::isnan(color_bounds[0]) || std::isnan(color_bounds[1])) {
    color_bounds = bounds_field(color_by)[0];
    color_bounds[1] += tol;
  }
  auto& elements = acc_mesh.elements();
  for (int i_elem = 0; i_elem < elements.size(); ++i_elem) {
    for (int i_contour = 0; i_contour < n_contour; ++i_contour) {
      double contour_val = contour_bounds[0] + (i_contour + 1)/(n_contour + 1.)*(contour_bounds[1] - contour_bounds[0]);
      auto& elem = elements[i_elem];
      // add contour line/surface
      otter_vis::add_contour(plt, elem, basis, contour, contour_val, n_div,
                             color_by, color_bounds, cmap_contour, transparent, status.flow_time, tol);
      // for 2d, color the flow field as well
      if (params.n_dim == 2) {
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
}
#endif

}
