#include <Boundary_condition.hpp>
#include <math.hpp>
#include <kernel_factory.hpp>
#include <Surface_rotation.hpp>

namespace cartdg
{

Freestream::Freestream(std::vector<double> freestream_state)
: fs{freestream_state}
{}

void Freestream::apply(Boundary_face& bf)
{
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  double* gf = bf.ghost_face();
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gf[i_var*nq + i_qpoint] = fs[i_var];
    }
  }
}

void Nonpenetration::apply(Boundary_face& bf)
{
  // fetch data
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  double* gh_f = bf.ghost_face();
  double* in_f = bf.inside_face();
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  // rotate into surface-based coordinates
  int i_dim = bf.i_dim();
  auto surf_rot = cartdg::kernel_factory<Surface_rotation>(params.n_dim, params.row_size, bf.jacobian_mat(), i_dim);
  surf_rot->to_surface(gh_f);
  // reflect normal component of momentum
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
    gh_f[i_dim*nq + i_qpoint] *= -1;
  }
  // rotate back into universal coordinates
  surf_rot->from_surface(gh_f);
}


void Copy::apply(Boundary_face& bf)
{
  double* in_f = bf.inside_face();
  double* gh_f = bf.ghost_face();
  for (int i_dof = 0; i_dof < bf.storage_params().n_dof()/bf.storage_params().row_size; ++i_dof) {
    gh_f[i_dof] = in_f[i_dof];
  }
}

void Nominal_pos::snap_vertices(Boundary_connection& con)
{
  const int stride = custom_math::pow(2, con.storage_params().n_dim - 1 - con.i_dim());
  for (int i_vert = 0; i_vert < con.storage_params().n_vertices(); ++i_vert) {
    if ((i_vert/stride)%2 == con.inside_face_sign()) {
      double pos = (con.element().nominal_position()[con.i_dim()] + con.inside_face_sign())*con.element().nominal_size();
      con.element().vertex(i_vert).pos[con.i_dim()] = pos;
    }
  }
}

void Surface_mbc::snap_vertices(Boundary_connection& con)
{
  const int stride = custom_math::pow(2, con.storage_params().n_dim - 1 - con.i_dim());
  for (int i_vert = 0; i_vert < con.storage_params().n_vertices(); ++i_vert) {
    if ((i_vert/stride)%2 == con.inside_face_sign()) {
      auto& pos = con.element().vertex(i_vert).pos;
      pos = sg->project_point(pos);
    }
  }
}

void Surface_mbc::snap_node_adj(Boundary_connection& con, const Basis& basis)
{
  if (!con.element().node_adjustments()) return; // Cartesian elements don't have `node_adjustments()`, so in this case just exit
  auto params {con.storage_params()};
  const int nfq = params.n_qpoint()/params.row_size;
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    // get position on this and opposite face
    std::array<double, 3> face_pos [2];
    for (int face_sign = 0; face_sign < 2; ++face_sign) {
      std::vector<double> pos {con.element().face_position(basis, 2*con.i_dim() + face_sign, i_qpoint)};
      for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) face_pos[face_sign][i_dim] = pos[i_dim];
      for (unsigned i_dim = pos.size(); i_dim < 3; ++i_dim) face_pos[face_sign][i_dim] = 0.;
    }
    // compute intersections
    auto sects = sg->line_intersections(face_pos[0], face_pos[1]);
    // set the node adjustment to match the nearest intersection, if any of them are reasonably close
    double best_adj = 0.;
    double distance = 1.;
    for (double sect : sects) {
      double adj = sect - con.inside_face_sign();
      if (std::abs(adj) < distance) {
        best_adj = adj;
        distance = std::abs(adj);
      }
    }
    con.element().node_adjustments()[(2*con.i_dim() + con.inside_face_sign())*nfq + i_qpoint] += best_adj;
  }
}

}
