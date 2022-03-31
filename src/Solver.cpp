#include <Solver.hpp>
#include <Mcs_cartesian.hpp>
#include <Mcs_deformed.hpp>
#include <Write_face.hpp>
#include <Neighbor_cartesian.hpp>
#include <Neighbor_deformed.hpp>
#include <Local_cartesian.hpp>
#include <Local_deformed.hpp>
#include <Tecplot_file.hpp>

namespace cartdg
{

void Solver::share_vertex_data(Element::shareable_value_access access_func, Vertex::reduction reduce)
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
  params{3, n_dim + 2, n_dim, row_size},
  acc_mesh{params, root_mesh_size},
  basis{row_size}
{}

void Solver::relax_vertices()
{
  // relax the vertices
  auto verts {acc_mesh.vertices()};
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].calc_relax();
  for (int i_vert = 0; i_vert < verts.size(); ++i_vert) verts[i_vert].apply_relax();
  // vertex relaxation will cause hanging vertices to drift away from the faces they are supposed to be coincident with
  // so now we put them back where they belong
  auto& matchers = acc_mesh.hanging_vertex_matchers();
  for (int i_match = 0; i_match < matchers.size(); ++i_match) {
    matchers[i_match].match(&Element::vertex_position<0>);
    matchers[i_match].match(&Element::vertex_position<1>);
    matchers[i_match].match(&Element::vertex_position<2>);
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
        if (dir.transpose()) { // already implies n_dim == 3
          Eigen::Map<Eigen::MatrixXd> face (elem_jac[1], rs, rs);
          face.transposeInPlace();
        }
        if (dir.flip_tangential()) {
          int free_dim = n_dim*(n_dim - 1)/2 - dir.i_dim[0] - dir.i_dim[1]; // which dimension is not involved in the connection?
          int stride = ((n_dim == 3) && (dir.i_dim[1] < free_dim)) ? rs : 1; // along what stride do we need to reverse the elements?
          for (int i_row = 0; i_row < nfq/rs; ++i_row) {
            Eigen::Map<Eigen::VectorXd, 0, Eigen::InnerStride<>> row (elem_jac[1] + i_row*rs/stride, rs, Eigen::InnerStride<>(stride));
            row.reverseInPlace();
          }
        }
        for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
          // take average of element face jacobians with appropriate axis permutations
          int col = j_dim;
          int tangential_sign = 1;
          // swap dir.i_dim[0] and dir.i_dim[1]
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
    Eigen::Map<Eigen::VectorXd> vert_tss (elem.vertex_time_step_scale(), params.n_vertices());
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
}

void Solver::update(double stability_ratio)
{
  const int nd = params.n_dim;
  const int rs = params.row_size;
  auto& elems = acc_mesh.elements();
  // compute time step
  double mcs = std::max((*kernel_factory<Mcs_cartesian>(nd, rs))(acc_mesh.cartesian().elements()),
                        (*kernel_factory<Mcs_deformed >(nd, rs))(acc_mesh.deformed ().elements()));
  double dt = stability_ratio*basis.max_cfl_convective()/params.n_dim/mcs;
  // record reference state for Runge-Kutta scheme
  const int n_dof = params.n_dof();
  for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
    double* state = elems[i_elem].stage(0);
    for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] = state[i_dof];
  }
  // perform update for each Runge-Kutta stage
  for (double rk_weight : {rk_weights[0]}) { // FIXME
    (*kernel_factory<Write_face>(nd, rs, basis))(elems);
    auto& bcs {acc_mesh.boundary_conditions()};
    auto& bc_cons {acc_mesh.boundary_connections()};
    for (int i_con = 0; i_con < bc_cons.size(); ++i_con) {
      for (int i_bc = 0; i_bc < bcs.size(); ++i_bc) {
        if (&bcs[i_bc] == bc_cons[i_con].boundary_condition()) bcs[i_bc].apply(bc_cons[i_con]);
      }
    }
    (*kernel_factory<Neighbor_cartesian>(nd, rs))(acc_mesh.cartesian().face_connections());
    (*kernel_factory<Neighbor_deformed >(nd, rs))(acc_mesh.deformed ().face_connections());
    (*kernel_factory<Local_cartesian>(nd, rs, basis, dt, rk_weight))(acc_mesh.cartesian().elements());
    (*kernel_factory<Local_deformed >(nd, rs, basis, dt, rk_weight))(acc_mesh.deformed ().elements());
  }
  // update status for reporting
  status.time_step = dt;
  status.flow_time += dt;
  ++status.iteration;
}

Iteration_status Solver::iteration_status()
{
  return status;
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
  std::vector<double> integral (integrand.n_var(params.n_dim));
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

void Solver::visualize_field(const Qpoint_func& output_variables, std::string name)
{
  const int n_dim = params.n_dim;
  const int n_vis = output_variables.n_var(n_dim); // number of variables to visualize
  const int n_sample = 20;
  const int n_qpoint = params.n_qpoint();
  const int n_corners {custom_math::pow(2, n_dim - 1)};
  const int nfqpoint = n_qpoint/params.row_size;
  const int n_block {custom_math::pow(n_sample, n_dim)};
  Eigen::MatrixXd interp {basis.interpolate(Eigen::VectorXd::LinSpaced(n_sample, 0., 1.))};
  std::vector<std::string> var_names;
  for (int i_vis = 0; i_vis < n_vis; ++i_vis) var_names.push_back(output_variables.variable_name(i_vis));
  Tecplot_file file {name, n_dim, var_names, status.flow_time};

  for (int i_elem = 0; i_elem < acc_mesh.elements().size(); ++i_elem)
  {
    Element& elem {acc_mesh.elements()[i_elem]};
    std::vector<double> pos (n_qpoint*n_dim);
    std::vector<double> to_vis (n_qpoint*n_vis);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      auto qpoint_vis = output_variables(elem, basis, i_qpoint, status.flow_time);
      for (int i_vis = 0; i_vis < n_vis; ++i_vis) {
        to_vis[i_vis*n_qpoint + i_qpoint] = qpoint_vis[i_vis];
      }
      auto qpoint_pos = elem.position(basis, i_qpoint);
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        pos[i_dim*n_qpoint + i_qpoint] = qpoint_pos[i_dim];
      }
    }
    // note: each visualization stage is enclosed in `{}` to ensure that only one
    // `Tecplot_file::Zone` is alive at a time

    // visualize edges
    if (n_dim > 1) // 1D elements don't really have edges
    {
      Tecplot_file::Line_segments edges {file, n_dim*n_corners, n_sample, "edges"};
      Eigen::MatrixXd boundary {basis.boundary()};
      for (int i_dim = 0; i_dim < n_dim; ++i_dim)
      {
        const int stride {custom_math::pow(params.row_size, n_dim - 1 - i_dim)};
        const int n_outer {n_qpoint/stride/params.row_size};

        auto extract_edge = [=](double* data, int n)
        {
          Eigen::MatrixXd edge {n_sample, n_corners*n};
          for (int i = 0; i < n; ++i) {
            Eigen::MatrixXd edge_qpoints {basis.row_size, n_corners};
            for (int i_qpoint = 0; i_qpoint < basis.row_size; ++i_qpoint) {
              Eigen::VectorXd qpoint_slab {nfqpoint};
              for (int i_outer = 0; i_outer < n_outer; ++i_outer) {
                for (int i_inner = 0; i_inner < stride; ++i_inner) {
                  qpoint_slab[i_outer*stride + i_inner] = data[i*n_qpoint + i_qpoint*stride + i_outer*stride*basis.row_size + i_inner];
                }
              }
              edge_qpoints.row(i_qpoint) = custom_math::hypercube_matvec(boundary, qpoint_slab);
            }
            for (int i_corner = 0; i_corner < n_corners; ++i_corner) {
              edge.col(i_corner*n + i) = interp*edge_qpoints.col(i_corner);
            }
          }
          return edge;
        };

        Eigen::MatrixXd edge_pos {extract_edge(pos.data(), n_dim)};
        Eigen::MatrixXd edge_state {extract_edge(to_vis.data(), n_vis)};
        for (int i_corner = 0; i_corner < n_corners; ++i_corner) {
          edges.write(edge_pos.data() + i_corner*n_dim*n_sample, edge_state.data() + i_corner*n_vis*n_sample);
        }
      }
    }

    { // visualize quadrature points
      Tecplot_file::Structured_block qpoints {file, basis.row_size, "element_qpoints"};
      qpoints.write(pos.data(), to_vis.data());
    }

    { // visualize interior (that is, quadrature point data interpolated to a fine mesh of sample points)
      Tecplot_file::Structured_block interior {file, n_sample, "element_interior"};
      Eigen::VectorXd interp_pos {n_block*n_dim};
      for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
        Eigen::Map<Eigen::VectorXd> qpoint_pos (pos.data() + i_dim*n_qpoint, n_qpoint);
        interp_pos.segment(i_dim*n_block, n_block) = custom_math::hypercube_matvec(interp, qpoint_pos);
      }
      Eigen::VectorXd interp_state {n_block*n_vis};
      for (int i_var = 0; i_var < n_vis; ++i_var) {
        Eigen::Map<Eigen::VectorXd> var (to_vis.data() + i_var*n_qpoint, n_qpoint);
        interp_state.segment(i_var*n_block, n_block) = custom_math::hypercube_matvec(interp, var);
      }
      interior.write(interp_pos.data(), interp_state.data());
    }
  }
}

}
