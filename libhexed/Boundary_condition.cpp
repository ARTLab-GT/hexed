#include <hexed/Boundary_condition.hpp>
#include <hexed/math.hpp>
#include <hexed/kernel_factory.hpp>
#include <hexed/constants.hpp>
#include <hexed/pde.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include <hexed/Printer.hpp>

namespace hexed {

void copy_state(Boundary_connection& con) {
  con.ghost().full_state() = con.inside().full_state();
}

void Flow_bc::apply_advection(Boundary_connection& con) {
  int nd = con.ghost().storage_params().n_dim;
  Array<double> inside_state = con.inside().advection_state();
  Array<double> ghost_state = con.ghost().advection_state();
  ghost_state(0, nd) = inside_state(0, nd);
  ghost_state(nd, end) = 2. - inside_state(nd, end);
}

void Flow_bc::apply_diffusion(Boundary_connection& con) {
  con.ghost().flow_state()(0) = con.inside().flow_state()(0);
}

void Flow_bc::flux_diffusion(Boundary_connection& con) {
  con.ghost().flow_state()(1) = -con.inside().flow_state()(1);
}

void Flow_bc::init_cache(Boundary_connection& con) {
  con.state_cache() = con.inside().flow_state()(0);
  con.flux_cache() = con.inside().flow_state()(1);
}

Freestream::Freestream(Mat<> freestream_state)
: fs{freestream_state}
{}

void Freestream::apply_state(Boundary_connection& con) {
  Array<double> ghost_state = con.ghost().flow_state()(0);
  auto params = con.ghost().storage_params();
  for (int i_var = 0; i_var < params.n_var; ++i_var) ghost_state(i_var) = fs(i_var);
}

Riemann_invariants::Riemann_invariants(Mat<> freestream_state)
: fs{freestream_state}
{}

template <int n_dim>
Mat<> apply_char(Mat<> state, Mat<> normal, int sign, Mat<> inside, Mat<> outside) {
  // compute characteristics
  typename pde::Navier_stokes<>::Pde<n_dim, 2>::Characteristics ch(state, normal);
  auto eigvals = ch.eigvals();
  auto decomp = ch.decomp(inside);
  auto fs_decomp = ch.decomp(outside);
  // set eigenvectors to inside or outside values depending on sign of eigenvalues
  for (int i_eig = 0; i_eig < 3; ++i_eig) {
    if (sign*eigvals(i_eig) > 0) decomp(Eigen::all, i_eig) = fs_decomp(Eigen::all, i_eig);
  }
  return decomp.rowwise().sum();
}

void Riemann_invariants::apply_state(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_state = con.ghost().flow_state()(0);
  Array<double> inside_state = con.inside().flow_state()(0);
  Array<double> normal = con.normal();
  int sign = 1 - 2*con.inside().sign(); // sign of velocity of incoming characteristics
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> state(params.n_var);
    // fetch data
    Mat<> inside(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) inside(i_var) = inside_state(i_var)[i_qpoint];
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = normal(i_dim)[i_qpoint];
    // compute characteristics
    switch (params.n_dim) {
      case 1:
        // set incoming characteristics to zero and leave outgoing alone
        state = apply_char<1>(inside, n, sign, inside, fs);
        break;
      case 2:
        state = apply_char<2>(inside, n, sign, inside, fs);
        break;
      case 3:
        state = apply_char<3>(inside, n, sign, inside, fs);
        break;
      default:
        throw std::runtime_error("invalid dimensionality");
    }
    // limit state to ensure thermodynamic admissibility
    state(params.n_dim) = std::max(state(params.n_dim), inside(params.n_dim)/2);
    double kin_ener = .5*state(Eigen::seqN(0, params.n_dim)).squaredNorm()/state(params.n_dim);
    double inside_kin_ener = .5*inside(Eigen::seqN(0, params.n_dim)).squaredNorm()/inside(params.n_dim);
    state(params.n_dim + 1) = std::max(kin_ener + std::max(state(params.n_dim + 1) - kin_ener,
                                                           (inside(params.n_dim + 1) - inside_kin_ener)/2), 0.);
    // write to ghost state
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      ghost_state(i_var)[i_qpoint] = state(i_var);
    }
  }
  // prime state cache with inside state
  con.state_cache() = inside_state();
}

void Riemann_invariants::apply_flux(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_flux = con.ghost().flow_state()(1);
  Array<double> inside_flux = con.inside().flow_state()(1);
  Array<double> state_cache = con.state_cache();
  Array<double> normal = con.normal();
  int sign = 1 - 2*con.inside().sign();
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    // fetch data
    Mat<> inside(params.n_var); // flux
    Mat<> cache(params.n_var); // state
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      inside(i_var) = inside_flux(i_var)[i_qpoint];
      cache(i_var) = state_cache(i_var)[i_qpoint];
    }
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = normal(i_dim)[i_qpoint];
    // compute characteristics
    Mat<> state;
    switch (params.n_dim) {
      case 1:
        state = apply_char<1>(cache, n, sign, Mat<>::Zero(params.n_var), inside); // set outgoing characteristic flux to zero and leave incoming alone
        break;
      case 2:
        state = apply_char<2>(cache, n, sign, Mat<>::Zero(params.n_var), inside);
        break;
      case 3:
        state = apply_char<3>(cache, n, sign, Mat<>::Zero(params.n_var), inside);
        break;
      default:
        throw std::runtime_error("invalid dimensionality");
    }
    // write to ghost flux
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      ghost_flux(i_var)[i_qpoint] = state(i_var);
    }
  }
}

void Pressure_outflow::apply_state(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_state = con.ghost().flow_state()(0);
  Array<double> inside_state = con.inside().flow_state()(0);
  Array<double> normal = con.normal();
  int sign = 2*con.inside().sign() - 1;
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    // fetch data
    Mat<> inside(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) inside(i_var) = inside_state(i_var)[i_qpoint];
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = normal(i_dim)[i_qpoint];
    Mat<> ghost = inside;
    // if subsonic, set pressure to specified value
    Mat<> mmtm = inside(Eigen::seqN(0, params.n_dim));
    double nrml_veloc = mmtm.dot(n)/inside(params.n_dim)/n.norm();
    double kin_ener = .5*mmtm.squaredNorm()/inside(params.n_dim);
    double pres = std::max(.4*(inside(params.n_dim + 1) - kin_ener), 0.);
    #if 1
    double sound_speed = std::sqrt(1.4*pres/inside(params.n_dim));
    if (nrml_veloc*sign < sound_speed) ghost(params.n_dim + 1) = pres_spec/.4 + kin_ener;
    #else
    double sound_speed = std::sqrt(std::abs(1.4*pres/inside(params.n_dim)));
    double nrml_mach = nrml_veloc*sign/sound_speed;
    double ramp_size = 0.1;
    double interp = std::min(0., std::max(1., (1. + ramp_size - nrml_mach)/ramp_size));
    ghost(params.n_dim + 1) += interp*(pres_spec/.4 + kin_ener - ghost(params.n_dim + 1));
    #endif
    // write to ghost flux
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      ghost_state(i_var)[i_qpoint] = ghost(i_var);
    }
  }
}

//! \todo make this formally well-posed
void Pressure_outflow::apply_flux(Boundary_connection& con) {
  con.ghost().flow_state()(1) = -con.inside().flow_state()(1);
}

void Freestream::apply_flux(Boundary_connection& con) {copy_state(con);}

void reflect_normal(double* gh_f, double* nrml, int nq, int nd) {
  for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
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

void reflect_momentum(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  Array<double> ghost_state = con.ghost().flow_state()(0);
  Array<double> inside_state = con.inside().flow_state()(0);
  ghost_state = inside_state;
  reflect_normal(ghost_state.data(), con.normal().data(), params.n_qpoint()/params.row_size, params.n_dim);
}

void Nonpenetration::apply_state(Boundary_connection& con) {
  reflect_momentum(con);
}

void Nonpenetration::apply_flux(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  Array<double> ghost_state = con.ghost().flow_state()(1);
  Array<double> inside_state = con.inside().flow_state()(1);
  ghost_state = -inside_state;
  reflect_normal(ghost_state.data(), con.normal().data(), params.n_qpoint()/params.row_size, params.n_dim);
}

void Nonpenetration::apply_advection(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_state = con.ghost().advection_state();
  Array<double> inside_state = con.inside().advection_state();
  Array<double> normal = con.normal();
  int nd = params.n_dim;
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> veloc(nd);
    Mat<> n(nd);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      n(i_dim) = normal(i_dim)[i_qpoint];
      veloc(i_dim) = inside_state(i_dim)[i_qpoint];
    }
    veloc -= 2*veloc.dot(n)*n/n.squaredNorm();
    for (int i_dim = 0; i_dim < nd; ++i_dim) ghost_state(i_dim)[i_qpoint] = veloc(i_dim);
  }
  ghost_state(nd, end) = inside_state(nd, end);
}

No_slip::No_slip(std::shared_ptr<Thermal_bc> thermal, double heat_rat, Transport_model visc,
                 Turbulence_model turb, double coercion)
: _coercion{coercion}
, _thermal{thermal}
, _viscosity{visc}
, _turb{turb}
, _heat_rat{heat_rat}
{}

double Thermal_equilibrium::ghost_heat_flux(Mat<> state, double) {
  double temp = state(last)*(heat_rat - 1)/state(state.size() - 2)/constants::specific_gas_air;
  double radiative_flux = emissivity*constants::stefan_boltzmann*math::pow(temp, 4);
  double conductive_flux = heat_transfer_coef*(temp - temperature);
  return radiative_flux + conductive_flux;
}

void No_slip::apply_state(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_state = con.ghost().flow_state()(0);
  Array<double> inside_state = con.inside().flow_state()(0);
  Array<double> state_cache = con.state_cache();
  Array<double> presc = con.prescribed_data();
  int nd = params.n_dim;
  // set ghost state
  // momentum
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    ghost_state(i_dim) = 2.*presc(i_dim)*inside_state(nd) - inside_state(i_dim);
  }
  // density
  ghost_state(nd) = inside_state(nd);
  // energy
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> state(params.n_dim + 2); // yes, this should be ignoring turbulence variables
    for (int i_var = 0; i_var < params.n_dim + 2; ++i_var) state(i_var) = inside_state(i_var)[i_qpoint];
    double kin_ener = .5*state(Eigen::seqN(0, nd)).squaredNorm()/state(nd);
    double internal_energy = state(last) - kin_ener;
    double ghost_energy = _thermal->ghost_energy(state);
    if (ghost_energy > internal_energy) ghost_energy += ghost_energy - internal_energy;
    else ghost_energy *= ghost_energy/internal_energy;
    ghost_state(nd + 1)[i_qpoint] = ghost_energy + kin_ener;
  }
  if (_turb == k_omega) {
    ghost_state(nd + 2) = -inside_state(nd + 2); // set turbulent kinetic energy to 0
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      // set dissipation based on wall roughness
      double mass = inside_state(nd)[i_qpoint];
      double energy = inside_state(nd + 1)[i_qpoint]/mass;
      for (int i_dim = 0; i_dim < nd; ++i_dim) {
        energy -= .5*math::pow(inside_state(i_dim)[i_qpoint]/mass, 2);
      }
      double dyn_visc = _viscosity.coefficient(std::sqrt(std::abs(energy*(_heat_rat - 1.)
                                                         /constants::specific_gas_air)));
      double roughness = presc(nd)[i_qpoint];
      double omega_wall = 4e4*dyn_visc/(mass*roughness*roughness);
      ghost_state(nd + 3)[i_qpoint] = 2*std::log(omega_wall)*mass - inside_state(params.n_dim + 3)[i_qpoint];
    }
  }
  // prime `state_cache` with average state for use in emissivity BC
  state_cache = (ghost_state + inside_state)/2.;
}

void No_slip::apply_flux(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  Array<double> ghost_state = con.ghost().flow_state()(1);
  Array<double> inside_state = con.inside().flow_state()(1);
  Array<double> normal = con.normal();
  Array<double> state_cache = con.state_cache();
  Array<double> presc = con.prescribed_data();
  // set momentum and mass flux (pretty straightforward)
  ghost_state(0, params.n_dim) = inside_state(0, params.n_dim);
  ghost_state(params.n_dim) = -inside_state(params.n_dim);
  // set energy flux depending on thermal boundary condition
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    double nrml = 0;
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      nrml += normal(i_dim)[i_qpoint]*normal(i_dim)[i_qpoint];
    }
    nrml = std::sqrt(nrml);
    int flux_sign = 2*con.inside().sign() - 1;
    Mat<> state(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) state(i_var) = state_cache(i_var)[i_qpoint];
    double inside_ener = inside_state(params.n_dim + 1)[i_qpoint];
    double ghost_heat = _thermal->ghost_heat_flux(state, inside_ener*flux_sign/nrml);
    ghost_state(params.n_dim + 1)[i_qpoint] = _coercion*(nrml*flux_sign*ghost_heat - inside_ener) + inside_ener;
  }
  // set turbulence variables
  ghost_state(params.n_dim + 2, end) = inside_state(params.n_dim + 2, end);
}

void No_slip::apply_advection(Boundary_connection& con) {
  auto params = con.ghost().storage_params();
  Array<double> ghost_state = con.ghost().advection_state();
  Array<double> inside_state = con.inside().advection_state();
  // set velocity to 0
  ghost_state(0, params.n_dim) = -inside_state(0, params.n_dim);
  // don't change advected scalar
  ghost_state(params.n_dim, end) = inside_state(params.n_dim, end);
}

void No_slip::set_prescribed(Interpreter& inter, Boundary_connection& con) {
  auto sub = inter.make_sub();
  auto params = con.ghost().storage_params();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint()/params.row_size;
  sub.variables->assign("pos", con.position());
  sub.variables->assign("normal", con.normal());
  auto expr = inter.variables->lookup<std::string>("wall_velocity");
  HEXED_ASSERT(expr, "`wall_velocity` must be specified for `No_slip`");
  sub.exec(expr.value());
  Array<double> data = con.prescribed_data();
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    sub.variables->assign_array(data(i_dim), "velocity" + std::to_string(i_dim));
  }
  if (_turb == k_omega) {
    Array<double> state = con.state_cache();
    Array<double> flux = con.flux_cache();
    Array<double> normal = con.normal();
    double area = con.inside().nominal_area();
    bool local = inter.variables->get<int>("local_roughness");
    if (local) {
      double max_rough = inter.variables->get<double>("hexed_max_roughness");
      double max_plus = inter.variables->get<double>("max_roughness_plus");
      for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
        double nrml = 0;
        double stress = 0;
        for (int i_dim = 0; i_dim < nd; ++i_dim) {
          stress += math::pow(flux(i_dim)[i_qpoint], 2);
          nrml += math::pow(normal(i_dim)[i_qpoint], 2);
        }
        stress = std::sqrt(stress/nrml)/area;
        double density = state(nd)[i_qpoint];
        double temperature = state(nd + 1)[i_qpoint]*(_heat_rat - 1)/density/constants::specific_gas_air;
        double dyn_visc = _viscosity.coefficient(std::sqrt(temperature));
        double friction_veloc = std::sqrt(stress/density);
        double inv_length = (friction_veloc*density)/dyn_visc;
        data(nd)[i_qpoint] = 1./std::pow(1./math::pow(max_rough, 4) + math::pow(inv_length/max_plus, 4), .25);
      }
    } else {
      data(nd) = inter.variables->get<double>("hexed_surface_roughness");
    }
  }
}

void Copy::apply_state(Boundary_connection& con) {
  copy_state(con);
}

void Copy::apply_flux(Boundary_connection& con) {
  copy_state(con);
}

void Copy::apply_advection(Boundary_connection& con) {
  copy_state(con);
}

void Outflow::apply_state(Boundary_connection& con) {
  copy_state(con);
}

void Outflow::apply_flux(Boundary_connection& con) {
  // set to negative of inside flux
  con.ghost().flow_state()(1) = -con.inside().flow_state()(1);
}

Expression_bc::Expression_bc(Interpreter& inter, std::string state_expr, std::string flux_expr)
: _inter{inter}, _exprs{state_expr, flux_expr}
{}

void Expression_bc::_apply(Boundary_connection& con, bool is_flux) {
  auto sub = _inter.make_sub();
  auto params = con.ghost().storage_params();
  for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
    sub.variables->assign("pos" + to_string(i_dim), con.position()(i_dim));
    sub.variables->assign("normal" + to_string(i_dim), con.normal()(i_dim));
  }
  Array<double> inside = con.inside().flow_state();
  std::array<std::string, 2> names {"state", "flux"};
  for (int i = 0; i < 2; ++i) {
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      sub.variables->assign(names[i] + to_string(i_var), inside(i)(i_var));
    }
  }
  sub.exec(_exprs[is_flux]);
  Array<double> ghost = con.ghost().flow_state()(is_flux);
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    sub.variables->assign_array(ghost(i_var), "ghost_" + names[is_flux] + to_string(i_var));
  }
}

}
