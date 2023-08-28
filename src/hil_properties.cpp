#include <hil_properties.hpp>

namespace hexed::hil_properties
{

void element(Namespace& space, Element& elem)
{
  space.assign("is_extruded", int(!elem.tree));
  space.assign("ref_level", elem.refinement_level());
  space.assign("nom_sz", elem.nominal_size());
  space.assign("res_bad", elem.resolution_badness);
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
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    space.assign("momentum" + std::to_string(i_dim), elem.stage(0)[i_dim*nq + i_qpoint]);
  }
  for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
    space.assign("momentum" + std::to_string(i_dim), 0.);
  }
  space.assign("mass", elem.stage(0)[params.n_dim*nq + i_qpoint]);
  space.assign("energy", elem.stage(0)[(params.n_dim + 1)*nq + i_qpoint]);
  space.assign("art_visc", elem.art_visc_coef()[i_qpoint]);
  space.assign("tss", elem.time_step_scale()[i_qpoint]);
}

}
