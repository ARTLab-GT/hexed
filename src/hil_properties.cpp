#include <hil_properties.hpp>

namespace hexed::hil_properties
{

void element(Namespace& space, Element& elem)
{
  space.assign("is_extruded", int(!elem.tree));
  space.assign("ref_level", elem.refinement_level());
  space.assign("nom_sz", elem.nominal_size());
  space.assign("uncertainty", elem.uncertainty);
  space.assign("snapping_problem", int(elem.snapping_problem));
  auto params = elem.storage_params();
  Eigen::Vector3d center;
  center.setZero();
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    center += elem.vertex(i_vert).pos;
  }
  center /= params.n_vertices();
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    space.assign("center" + std::to_string(i_dim), center(i_dim));
  }
}

void position(Namespace& space, Element& elem, const Basis& basis, int i_qpoint)
{
  auto pos = elem.position(basis, i_qpoint);
  for (unsigned i_dim = 0; i_dim < pos.size(); ++i_dim) {
    space.assign("pos" + std::to_string(i_dim), pos[i_dim]);
  }
  for (int i_dim = pos.size(); i_dim < 3; ++i_dim) {
    space.assign("pos" + std::to_string(i_dim), 0.);
  }
}

void state(Namespace& space, Element& elem, int i_qpoint)
{
  auto params = elem.storage_params();
  int nq = params.n_qpoint();
  std::string prefixes [2] {{""}, {"prev_"}};
  for (int i_stage = 0; i_stage < 2; ++i_stage) {
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      space.assign(prefixes[i_stage] + "momentum" + std::to_string(i_dim), elem.stage(i_stage)[i_dim*nq + i_qpoint]);
    }
    for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
      space.assign(prefixes[i_stage] + "momentum" + std::to_string(i_dim), 0.);
    }
    space.assign(prefixes[i_stage] + "density", elem.stage(i_stage)[params.n_dim*nq + i_qpoint]);
    space.assign(prefixes[i_stage] + "energy", elem.stage(i_stage)[(params.n_dim + 1)*nq + i_qpoint]);
  }
  space.assign("art_visc", elem.art_visc_coef()[i_qpoint]);
  space.assign("tss", elem.time_step_scale()[i_qpoint]);
}

void surface(Namespace& space, Boundary_connection& con, int i_fqpoint)
{
  bool viscous = space.lookup<std::string>("viscosity_model").value() != "" || space.lookup<std::string>("conductivity_model").value() != "";;
  auto params = con.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  // fetch surface normal
  int nrml_sign = 1 - 2*con.inside_face_sign();
  Mat<> nrml(params.n_dim);
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    nrml(i_dim) = con.surface_normal()[i_dim*nfq + i_fqpoint];
  }
  double nrml_mag = nrml.norm();
  std::vector<double> flux;
  std::vector<double> state;
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    // fetch momentum flux in reference space
    double ref_flux = -con.flux_cache()[i_var*nfq + i_fqpoint];
    // compute stress
    flux.push_back((nrml_mag > 1e-3 && viscous) ? ref_flux*nrml_sign/nrml_mag : 0.);
    state.push_back(con.inside_face(false)[i_var*nfq + i_fqpoint]);
  }
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    space.assign("normal" + std::to_string(i_dim), nrml(i_dim)*nrml_sign/nrml_mag);
    space.assign("visc_stress" + std::to_string(i_dim), flux[i_dim]);
    space.assign("momentum" + std::to_string(i_dim), state[i_dim]);
  }
  for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
    space.assign("normal" + std::to_string(i_dim), 0.);
    space.assign("visc_stress" + std::to_string(i_dim), 0.);
    space.assign("momentum" + std::to_string(i_dim), 0.);
  }
  space.assign("density", state[params.n_dim]);
  space.assign("energy", state[params.n_dim + 1]);
  space.assign("mass_flux", flux[params.n_dim]);
  space.assign("heat_flux", flux[params.n_dim + 1]);
}

}
