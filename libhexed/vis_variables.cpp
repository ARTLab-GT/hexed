#include <hexed/vis_variables.hpp>
#include <hexed/Tree.hpp>
#include <Gauss_legendre.hpp>

namespace hexed::vis_variables {

std::string index(std::string name, int i) {return name + std::to_string(i);}

void element(Namespace& space, Element& elem) {
  space.assign("is_extruded", int(elem.is_extruded()));
  space.assign("ref_level", elem.refinement_level());
  space.assign("aniso_ref_level", elem.aniso_ref_level());
  space.assign("mask", elem.mask());
  space.assign("nom_sz", elem.nominal_size());
  space.assign("wall_distance", elem.wall_distance());
  space.assign("wall_dimension", elem.wall_dimension());
  space.assign("has_wall", int(elem.has_wall()));
  space.assign("uncertainty", elem.uncertainty);
  space.assign("snapping_problem", int(elem.snapping_problem));
  space.assign("is_deformed", int(elem.deformed()));
  space.assign("rms_residual", elem.residual);
  auto params = elem.storage_params();
  Mat<3> center;
  center.setZero();
  int sharp = 0;
  for (int i_vert = 0; i_vert < params.n_vertices(); ++i_vert) {
    center += elem.shape().vertex(i_vert).point({});
    next::Vertex& v = elem.active_shape().vertex(i_vert);
    sharp = sharp || v.snapped_edge >= 0 || v.snapped_point >= 0;
  }
  space.assign("sharp", sharp);
  center /= params.n_vertices();
  for (int i_dim = 0; i_dim < 3; ++i_dim) {
    space.assign(index("center", i_dim), center(i_dim));
    space.assign(index("aniso_ref_level", i_dim),
                 i_dim < params.n_dim ? elem.tree.value().anisotropic_refinement_level()[i_dim] : 0);
    space.assign(index("nominal_shape", i_dim), i_dim < params.n_dim ? elem.nominal_shape(i_dim) : 0.);
    space.assign(index("spectral_uncertainty", i_dim), i_dim < params.n_dim ? elem.spectral_uncert()[i_dim] : 0);
    space.assign("flux_uncertainty", elem.flux_uncert);
    space.assign<int>(index("sharp", i_dim), elem.is_sharp(i_dim));
  }
  space.assign("max_bulk_art_visc", Array<double>({params.n_qpoint()}, elem.bulk_av_coef()).extreme(1));
  Array<double> flow_state = elem.flow_state().copy();
  Array<double> spec_int_ener = flow_state(params.n_dim + 1)/flow_state(params.n_dim);
  Array<double> velocity = flow_state(0, params.n_dim);
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    velocity(i_dim) /= flow_state(params.n_dim);
    spec_int_ener -= .5*velocity(i_dim)*velocity(i_dim);
  }
  double sound_speed = std::sqrt(1.4*0.4*spec_int_ener.extreme(0));
  Gauss_legendre basis(params.row_size);
  Array<double> position = elem.position(basis);
  Mat<dyn, dyn> boundary = basis.boundary();
  int sonic = 0;
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    Array<double> face_veloc({params.n_dim, 2, params.n_face_qpoint()});
    Array<double> face_pos({params.n_dim, 2, params.n_face_qpoint()});
    for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) {
      Mat<> veloc = velocity(j_dim).vector();
      Mat<> pos = position(j_dim).vector();
      for (int sign = 0; sign < 2; ++sign) {
        face_veloc(j_dim)(sign).vector() = math::dimension_matvec(boundary(sign, all), veloc, i_dim);
        face_pos(j_dim)(sign).vector() = math::dimension_matvec(boundary(sign, all), pos, i_dim);
      }
    }
    for (int i_qpoint = 0; i_qpoint < params.n_face_qpoint(); ++i_qpoint) {
      for (int j_qpoint = 0; j_qpoint < params.n_face_qpoint(); ++j_qpoint) {
        int k_qpoint = i_dim == elem.wall_dimension() ? i_qpoint : j_qpoint;
        Mat<> direction(params.n_dim);
        for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) {
          direction(j_dim) = face_pos(j_dim)(1)[k_qpoint] - face_pos(j_dim)(0)[i_qpoint];
        }
        direction.normalize();
        double veloc0 = 0;
        double veloc1 = 0;
        for (int j_dim = 0; j_dim < params.n_dim; ++j_dim) {
          veloc0 += face_veloc(j_dim)(0)[i_qpoint]*direction(j_dim);
          veloc1 += face_veloc(j_dim)(1)[k_qpoint]*direction(j_dim);
        }
        sonic = sonic || (veloc0 >  sound_speed*1.03 && veloc1 <  sound_speed*0.97);
        sonic = sonic || (veloc0 > -sound_speed*0.97 && veloc1 < -sound_speed*1.03);
      }
    }
  }
  space.assign<int>("sonic", sonic);
  space.assign<int>("has_shock", elem.has_shock);
  space.assign<int>("time_step_ratio_diffusion", elem.ts_ratio_diffusion);
  space.assign<int>("time_step_ratio_decay", elem.ts_ratio_decay);
  space.assign<double>("local_av_width", Array<double>({params.n_qpoint()}, elem.laplacian_av_coef()).extreme(0));
}

void position(Namespace& space, Element& elem, const Basis& basis) {
  Array<double> pos {elem.position(basis)};
  for (int i_dim = 0; i_dim < pos.shape()[0]; ++i_dim) {
    space.assign(index("pos", i_dim), pos(i_dim).copy());
  }
  for (int i_dim = pos.shape()[0]; i_dim < 3; ++i_dim) {
    space.assign(index("pos", i_dim), 0.);
  }
  double* jac = elem.jacobian_determinant();
  if (jac) {
    space.assign("jacobian_det", Array<double>({elem.storage_params().n_qpoint()}, jac));
  } else {
    space.assign("jacobian_det", 1);
  }
}

void state(Namespace& space, Element& elem) {
  auto params = elem.storage_params();
  int nq = params.n_qpoint();
  auto assign_state = [&](std::string name, int i_var) {
    space.assign(name, Array<double>({nq}, elem.state() + i_var*nq).copy());
    space.assign("residual_" + name, Array<double>({nq}, elem.residual_cache() + i_var*nq).copy());
  };
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) assign_state("momentum" + std::to_string(i_dim), i_dim);
  for (int i_dim = params.n_dim; i_dim < 3; ++i_dim) {
    space.assign("momentum" + std::to_string(i_dim), 0.);
    space.assign("residual_momentum" + std::to_string(i_dim), 0.);
  }
  assign_state("density", params.n_dim);
  assign_state("energy", params.n_dim + 1);
  if (params.n_var == params.n_dim + 4) {
    assign_state("turbulent_kinetic_energy", params.n_dim + 2);
    assign_state("turbulent_dissipation_bassi", params.n_dim + 3);
  }
  space.assign("bulk_art_visc", Array<double>({nq}, elem.bulk_av_coef()));
  space.assign("laplacian_art_visc", Array<double>({nq}, elem.laplacian_av_coef()));
  space.assign("tss", Array<double>({nq}, elem.time_step_scale()));
  for (int i_var = 0; i_var < config::debug_variables; ++i_var) {
    Array<double> data({nq}, elem.debug_variables() + i_var*params.n_qpoint());
    space.assign("debug_var" + std::to_string(i_var), data());
  }
  if (space.get<int>("vis_art_visc_vars")) {
    int n_adv = params.n_offset*params.n_advection(params.row_size);
    Array<double> advection({n_adv, nq}, elem.advection_state());
    for (int i_adv = 0; i_adv < n_adv; ++i_adv) {
      space.assign("av_advection" + to_string(i_adv), advection(i_adv).copy());
    }
    Array<double> forcing({params.n_forcing, nq}, elem.art_visc_forcing());
    for (int i_force = 0; i_force < params.n_forcing; ++i_force) {
      space.assign("av_forcing" + to_string(i_force), forcing(i_force).copy());
    }
  }
}

void field(Namespace& space, Element& elem, const Basis& b) {
  element(space, elem);
  position(space, elem, b);
  state(space, elem);
}

void surface(Namespace& space, Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  // fetch surface data
  int nrml_sign = 1 - 2*con.inside().sign();
  // copying prevents these things from being inadvertently modified
  Array<double> nrml = con.normal().copy();
  Array<double> pos = con.position().copy();
  Array<double> state = con.inside().flow_state()(0).copy();
  Array<double> flux = con.inside().flow_state()(1).copy();
  Array<double> ref_flux = con.flux_cache().copy();
  double area = con.inside().nominal_area();
  // compute normals and fluxes
  for (int i_fqpoint = 0; i_fqpoint < nfq; ++i_fqpoint) {
    double norm = 0;
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) norm += math::pow(nrml(i_dim)[i_fqpoint], 2);
    norm = std::sqrt(norm);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) nrml(i_dim)[i_fqpoint] /= nrml_sign*norm;
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      flux(i_var)[i_fqpoint] = norm > 1e-6 ? -ref_flux(i_var)[i_fqpoint]*nrml_sign/norm/area : 0;
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
  if (params.n_var >= params.n_dim + 4) {
    space.assign("turbulent_kinetic_energy", state(params.n_dim + 2).copy());
    space.assign("turbulent_dissipation_bassi", state(params.n_dim + 3).copy());
    space.assign("roughness_height", con.prescribed_data()(params.n_dim).copy());
  } else {
    space.assign("roughness_height", 0.);
  }
  Element* elem = con.inside().element();
  HEXED_ASSERT(elem, "Boundary face has no element.")
  element(space, *elem);
  space.assign("wall_spacing", space.get<double>("wall_distance"));
}

}
