#include <hexed/vis_variables.hpp>

namespace hexed::vis_variables {

std::string index(std::string name, int i) {return name + std::to_string(i);}

void element(Namespace& space, Element& elem) {
  space.assign("is_extruded", int(!elem.tree));
  space.assign("ref_level", elem.refinement_level());
  space.assign("aniso_ref_level", elem.aniso_ref_level());
  space.assign("mask", elem.mask());
  space.assign("nom_sz", elem.nominal_size());
  space.assign("uncertainty", elem.uncertainty);
  space.assign("snapping_problem", int(elem.snapping_problem));
  auto params = elem.storage_params();
  Mat<3> center;
  center.setZero();
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    center += elem.shape().vertex(i_vert).point({});
  }
  center /= params.n_vertices();
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    space.assign(index("center", i_dim), center(i_dim));
  }
}

void position(Namespace& space, Element& elem, const Basis& basis) {
  Array<double> pos {elem.position(basis)};
  for (int i_dim = 0; i_dim < pos.shape()[0]; ++i_dim) {
    space.assign(index("pos", i_dim), pos(i_dim).copy());
  }
  for (int i_dim = pos.shape()[0]; i_dim < 3; ++i_dim) {
    space.assign(index("pos", i_dim), 0.);
  }
}

void state(Namespace& space, Element& elem) {
  auto params = elem.storage_params();
  int nq = params.n_qpoint();
  auto assign_state = [&](std::string name, int i_var) {
    space.assign(name, Array<double>({nq}, elem.state() + i_var*nq));
    space.assign("residual_" + name, Array<double>({nq}, elem.residual_cache() + i_var*nq));
  };
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) assign_state("momentum" + std::to_string(i_dim), i_dim);
  for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
    space.assign("momentum" + std::to_string(i_dim), 0.);
    space.assign("residual_momentum" + std::to_string(i_dim), 0.);
  }
  assign_state("density", params.n_dim);
  assign_state("energy", params.n_dim + 1);
  if (params.n_var == params.n_dim + 5) {
    assign_state("turbulent_kinetic_energy", params.n_dim + 2);
    assign_state("turbulent_dissipation_bassi", params.n_dim + 3);
    assign_state("turbulent_production", params.n_dim + 4);
  }
  space.assign("bulk_art_visc", Array<double>({nq}, elem.bulk_av_coef()));
  space.assign("laplacian_art_visc", Array<double>({nq}, elem.laplacian_av_coef()));
  space.assign("tss", Array<double>({nq}, elem.time_step_scale()));
  for (int i_var = 0; i_var < config::debug_variables; ++i_var) {
    Array<double> data({nq}, elem.debug_variables() + i_var*params.n_qpoint());
    space.assign("debug_var" + std::to_string(i_var), data());
  }
}

void field(Namespace& space, Element& elem, const Basis& b) {
  element(space, elem);
  position(space, elem, b);
  state(space, elem);
}

void surface(Namespace& space, Boundary_connection& con) {
  auto params = con.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  // fetch surface data
  int nrml_sign = 1 - 2*con.inside_face_sign();
  Array<double> nrml {Array<double>({params.n_dim, nfq}, con.surface_normal()).copy()};
  Array<double> pos({params.n_dim, nfq}, con.surface_position());
  Array<double> state({params.n_var, nfq}, con.inside_face(false));
  Array<double> flux({params.n_var, nfq});
  double* ref_flux = con.flux_cache();
  // compute normals and fluxes
  for (int i_fqpoint = 0; i_fqpoint < nfq; ++i_fqpoint) {
    double norm = 0;
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) norm += math::pow(nrml(i_dim)[i_fqpoint], 2);
    norm = std::sqrt(norm);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) nrml(i_dim)[i_fqpoint] /= nrml_sign*norm;
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      flux(i_var)[i_fqpoint] = norm > 1e-6 ? -ref_flux[i_var*nfq + i_fqpoint]*nrml_sign/norm : 0;
    }
  }
  // assign variables
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    space.assign(index("pos", i_dim), pos(i_dim).copy());
    space.assign(index("normal", i_dim), nrml(i_dim).copy());
    space.assign(index("visc_stress", i_dim), flux(i_dim).copy());
    space.assign(index("momentum", i_dim), state(i_dim).copy());
  }
  for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
    space.assign(index("pos", i_dim), 0.);
    space.assign(index("normal", i_dim), 0.);
    space.assign(index("visc_stress", i_dim), 0.);
    space.assign(index("momentum", i_dim), 0.);
  }
  space.assign("density", state(params.n_dim).copy());
  space.assign("energy", state(params.n_dim + 1).copy());
  space.assign("mass_flux", flux(params.n_dim).copy());
  space.assign("heat_flux", flux(params.n_dim + 1).copy());
}

}
