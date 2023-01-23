#include <Boundary_condition.hpp>
#include <math.hpp>
#include <kernel_factory.hpp>

namespace hexed
{

void Flow_bc::apply_advection(Boundary_face& bf)
{
  auto params = bf.storage_params();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint()/params.row_size;
  double* in_f = bf.inside_face();
  double* gh_f = bf.ghost_face();
  // set velocity equal to inside
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gh_f[i_dim*nq + i_qpoint] = in_f[i_dim*nq + i_qpoint];
    }
  }
  // set advected scalar to initial value
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
    gh_f[nd*nq + i_qpoint] = 1.;
  }
}

Freestream::Freestream(std::vector<double> freestream_state)
: fs{freestream_state}
{}

void copy_state(Boundary_face& bf)
{
  double* in_f = bf.inside_face();
  double* gh_f = bf.ghost_face();
  int n_face_dof = bf.storage_params().n_dof()/bf.storage_params().row_size;
  for (int i_dof = 0; i_dof < n_face_dof; ++i_dof) {
    gh_f[i_dof] = in_f[i_dof];
    gh_f[i_dof + 2*n_face_dof] = in_f[i_dof + 2*n_face_dof];
  }
}

void Freestream::apply_state(Boundary_face& bf)
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

Function_bc::Function_bc(const Surface_func& func_arg) : func{func_arg} {}

void Function_bc::apply_state(Boundary_face& bf)
{
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = params.n_var;
  double* gh_f = bf.ghost_face();
  double* in_f = bf.inside_face();
  double* nrml = bf.surface_normal();
  double* pos = bf.surface_position();
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
    // fetch shared/inside face data
    std::vector<double> n(nd);
    std::vector<double> p(nd);
    std::vector<double> s(nv);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      n[i_dim] = nrml[i_dim*nq + i_qpoint];
      p[i_dim] =  pos[i_dim*nq + i_qpoint];
    }
    for (int i_var = 0; i_var < nv; ++i_var) s[i_var] = in_f[i_var*nq + i_qpoint];
    // apply function
    auto state = func(p, 0, s, n);
    // write result to ghost face
    for (int i_var = 0; i_var < nv; ++i_var) {
      gh_f[i_var*nq + i_qpoint] = state[i_var];
    }
  }
}

void Function_bc::apply_flux(Boundary_face& bf)
{
  copy_state(bf);
}

void Freestream::apply_flux(Boundary_face& bf)
{
  copy_state(bf);
}

void reflect_normal(double* gh_f, double* nrml, int nq, int nd)
{
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint)
  {
    double dot = 0.;
    double norm_sq = 0.;
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      double n = nrml[i_dim*nq + i_qpoint];
      dot += gh_f[i_dim*nq + i_qpoint]*n;
      norm_sq += n*n;
    }
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      gh_f[i_dim*nq + i_qpoint] -= 2*dot*nrml[i_dim*nq + i_qpoint]/norm_sq;
    }
  }
}

void reflect_momentum(Boundary_face& bf)
{
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face();
  double* in_f = bf.inside_face();
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  reflect_normal(gh_f, bf.surface_normal(), params.n_qpoint()/params.row_size, params.n_dim);
}

void Nonpenetration::apply_state(Boundary_face& bf)
{
  reflect_momentum(bf);
}

void Nonpenetration::apply_flux(Boundary_face& bf)
{
  // fetch data
  auto params = bf.storage_params();
  // FIXME: the first iteration here is deprecated
  for (int offset : {0, 2*params.n_dof()/params.row_size}) {
    double* gh_f = bf.ghost_face() + offset;
    double* in_f = bf.inside_face() + offset;
    // initialize to negative of inside flux
    for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
    // un-invert normal component of momentum
    reflect_normal(gh_f, bf.surface_normal(), params.n_qpoint()/params.row_size, params.n_dim);
  }
}

void Nonpenetration::apply_advection(Boundary_face& bf)
{
  reflect_momentum(bf);
}

No_slip::No_slip(Thermal_type type, double value) : t{type}, v{value} {}

void No_slip::apply_state(Boundary_face& bf)
{
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face();
  double* in_f = bf.inside_face();
  int nfq = params.n_qpoint()/params.row_size;
  for (int i_dof = 0; i_dof < params.n_dim*nfq; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
  for (int i_dof = params.n_dim*nfq; i_dof < (params.n_dim + 1)*nfq; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  for (int i_dof = (params.n_dim + 1)*nfq; i_dof < (params.n_dim + 2)*nfq; ++i_dof) {
    gh_f[i_dof] = (t == internal_energy) ? 2*v - in_f[i_dof] : in_f[i_dof];
  }
}

void No_slip::apply_flux(Boundary_face& bf)
{
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  for (int offset : {0, 2*params.n_dof()/params.row_size}) {
    double* gh_f = bf.ghost_face() + offset;
    double* in_f = bf.inside_face() + offset;
    for (int i_dof = 0; i_dof < params.n_dim*nfq; ++i_dof) gh_f[i_dof] = in_f[i_dof];
    for (int i_dof = params.n_dim*nfq; i_dof < (params.n_dim + 1)*nfq; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
    for (int i_dof = (params.n_dim + 1)*nfq; i_dof < (params.n_dim + 2)*nfq; ++i_dof) {
      gh_f[i_dof] = (t == heat_flux) ? 2*v - in_f[i_dof] : in_f[i_dof];
    }
  }
}

void No_slip::apply_advection(Boundary_face& bf)
{
  apply_state(bf);
}

void Copy::apply_state(Boundary_face& bf)
{
  copy_state(bf);
}

void Copy::apply_flux(Boundary_face& bf)
{
  copy_state(bf);
}

void Copy::apply_advection(Boundary_face& bf)
{
  copy_state(bf);
}

void Outflow::apply_state(Boundary_face& bf)
{
  copy_state(bf);
}

void Outflow::apply_flux(Boundary_face& bf)
{
  // set to negative of inside flux
  auto params = bf.storage_params();
  // FIXME: the first iteration here is deprecated
  for (int offset : {0, 2*params.n_dof()/params.row_size}) {
    double* gh_f = bf.ghost_face() + offset;
    double* in_f = bf.inside_face() + offset;
    for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
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

Surface_set::Surface_set(std::vector<Surface_geometry*> ptrs)
{
  for (Surface_geometry* ptr : ptrs) geoms.emplace_back(ptr);
}

std::array<double, 3> Surface_set::project_point(std::array<double, 3> point)
{
  std::array<double, 3> proj;
  double dist_sq = std::numeric_limits<double>::max();
  for (auto& geom : geoms) {
    auto p = geom->project_point(point);
    double d = 0;
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      d += (p[i_dim] - point[i_dim])*(p[i_dim] - point[i_dim]);
    }
    if (d < dist_sq) {
      proj = p;
      dist_sq = d;
    }
  }
  return proj;
}

std::vector<double> Surface_set::line_intersections(std::array<double, 3> point0, std::array<double, 3> point1)
{
  std::vector<double> inter;
  for (auto& geom : geoms) {
    auto geom_inter = geom->line_intersections(point0, point1);
    inter.insert(inter.end(), geom_inter.begin(), geom_inter.end());
  }
  return inter;
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
