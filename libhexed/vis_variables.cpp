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
  Eigen::Vector3d center;
  center.setZero();
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    center += elem.vertex(i_vert).pos;
  }
  center /= params.n_vertices();
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    space.assign(index("center", i_dim), center(i_dim));
  }
}

void position(Namespace& space, Element& elem, const Basis& basis) {
  Array<double> pos {elem.vis_position(basis)};
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
  space.assign("bulk_art_visc", Array<double>({nq}, elem.bulk_av_coef()));
  space.assign("laplacian_art_visc", Array<double>({nq}, elem.laplacian_av_coef()));
  space.assign("tss", Array<double>({nq}, elem.time_step_scale()));
}

}
