#include <Boundary_condition.hpp>
#include <connection.hpp>
#include <math.hpp>
#include <kernel_factory.hpp>
#include <constants.hpp>
#include <pde.hpp>
#include <Gauss_lobatto.hpp>

namespace hexed {

void copy_state(Boundary_face& bf) {
  int n_face_dof = bf.storage_params().n_dof()/bf.storage_params().row_size;
  for (bool is_ldg : {0, 1}) {
    double* in_f = bf.inside_face(is_ldg);
    double* gh_f = bf.ghost_face(is_ldg);
    for (int i_dof = 0; i_dof < n_face_dof; ++i_dof) {
      gh_f[i_dof] = in_f[i_dof];
    }
  }
}

void Flow_bc::apply_advection(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint()/params.row_size;
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  // set velocity equal to inside
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gh_f[i_dim*nq + i_qpoint] = in_f[i_dim*nq + i_qpoint];
    }
  }
  // set advected scalar to initial value
  for (int i_adv = 0; i_adv < params.n_advection(params.row_size)*nq; ++i_adv) {
    gh_f[nd*nq + i_adv] = 2. - in_f[nd*nq + i_adv];
  }
}

void Flow_bc::apply_diffusion(Boundary_face& bf) {
  auto params = bf.storage_params();
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  // set state equal to inside
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) {
    gh_f[i_dof] = in_f[i_dof];
  }
}

void Flow_bc::flux_diffusion(Boundary_face& bf) {
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face(true);
  double* in_f = bf.inside_face(true);
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
}

Freestream::Freestream(Mat<> freestream_state)
: fs{freestream_state}
{}

void Freestream::apply_state(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  double* gf = bf.ghost_face(false);
  for (int i_var = 0; i_var < params.n_var; ++i_var) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gf[i_var*nq + i_qpoint] = fs(i_var);
    }
  }
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

void Riemann_invariants::apply_state(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  double* nrml = bf.surface_normal();
  int sign = 1 - 2*bf.inside_face_sign(); // sign of velocity of incoming characteristics
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> state;
    // fetch data
    Mat<> inside(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) inside(i_var) = in_f[i_var*nfq + i_qpoint];
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = nrml[i_dim*nfq + i_qpoint];
    // compute characteristics
    switch (params.n_dim) {
      case 1:
        state = apply_char<1>(inside, n, sign, inside, fs); // set incoming characteristics to zero and leave outgoing alone
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
    state(params.n_dim + 1) = std::max(kin_ener + std::max(state(params.n_dim + 1) - kin_ener, (inside(params.n_dim + 1) - inside_kin_ener)/2), 0.);
    // write to ghost state
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      gh_f[i_var*nfq + i_qpoint] = state(i_var);
    }
  }
  // prime state cache with inside state
  double* sc = bf.state_cache();
  for (int i_dof = 0; i_dof < params.n_var*nfq; ++i_dof) {
    sc[i_dof] = in_f[i_dof];
  }
}

void Riemann_invariants::apply_flux(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  double* sc = bf.state_cache();
  double* in_f = bf.inside_face(true);
  double* gh_f = bf.ghost_face(true);
  double* nrml = bf.surface_normal();
  int sign = 1 - 2*bf.inside_face_sign();
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    // fetch data
    Mat<> inside(params.n_var); // flux
    Mat<> cache(params.n_var); // state
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      inside(i_var) = in_f[i_var*nfq + i_qpoint];
      cache(i_var) = sc[i_var*nfq + i_qpoint];
    }
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = nrml[i_dim*nfq + i_qpoint];
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
      gh_f[i_var*nfq + i_qpoint] = state(i_var);
    }
  }
}

void Pressure_outflow::apply_state(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nfq = params.n_qpoint()/params.row_size;
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  double* nrml = bf.surface_normal();
  int sign = 2*bf.inside_face_sign() - 1;
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    // fetch data
    Mat<> inside(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      inside(i_var) = in_f[i_var*nfq + i_qpoint];
    }
    Mat<> n(params.n_dim);
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) n(i_dim) = nrml[i_dim*nfq + i_qpoint];
    Mat<> ghost = inside;
    // if subsonic, set pressure to specified value
    Mat<> mmtm = inside(Eigen::seqN(0, params.n_dim));
    double nrml_veloc = mmtm.dot(n)/inside(params.n_dim)/n.norm();
    double kin_ener = .5*mmtm.squaredNorm()/inside(params.n_dim);
    double pres = std::max(.4*(inside(params.n_dim + 1) - kin_ener), 0.);
    double sound_speed = std::sqrt(1.4*pres/inside(params.n_dim));
    if (nrml_veloc*sign < sound_speed) ghost(params.n_dim + 1) = pres_spec/.4 + kin_ener;
    // write to ghost flux
    for (int i_var = 0; i_var < params.n_var; ++i_var) {
      gh_f[i_var*nfq + i_qpoint] = ghost(i_var);
    }
  }
}

//! \todo make this formally well-posed
void Pressure_outflow::apply_flux(Boundary_face& bf) {
  // set to negative of inside flux
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face(true);
  double* in_f = bf.inside_face(true);
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
}

Function_bc::Function_bc(const Surface_func& func_arg) : func{func_arg} {}

void Function_bc::apply_state(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = params.n_var;
  double* gh_f = bf.ghost_face(false);
  double* in_f = bf.inside_face(false);
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

void Cache_bc::apply_state(Boundary_face& bf) {
  int n_face_dof = bf.storage_params().n_dof()/bf.storage_params().row_size;
  double* sc = bf.state_cache();
  double* gh_state = bf.ghost_face(false);
  double* in_flux = bf.inside_face(true);
  double* gh_flux = bf.ghost_face(true);
  for (int i_dof = 0; i_dof < n_face_dof; ++i_dof) {
    gh_state[i_dof] = sc[i_dof];
    gh_flux[i_dof] = in_flux[i_dof];
  }
}

void Cache_bc::init_cache(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nq = params.n_qpoint()/params.row_size;
  const int nd = params.n_dim;
  const int nv = params.n_var;
  double* cache = bf.state_cache();
  double* in_f = bf.inside_face(false);
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
    auto state = (*func)(p, 0, s, n);
    // write result to ghost face
    for (int i_var = 0; i_var < nv; ++i_var) {
      cache[i_var*nq + i_qpoint] = state[i_var];
    }
  }
}

void Function_bc::apply_flux(Boundary_face& bf) {copy_state(bf);}
void Freestream::apply_flux(Boundary_face& bf) {copy_state(bf);}
void Cache_bc::apply_flux(Boundary_face& bf) {copy_state(bf);}

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

void reflect_momentum(Boundary_face& bf) {
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face(false);
  double* in_f = bf.inside_face(false);
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  reflect_normal(gh_f, bf.surface_normal(), params.n_qpoint()/params.row_size, params.n_dim);
}

void Nonpenetration::apply_state(Boundary_face& bf) {
  reflect_momentum(bf);
}

void Nonpenetration::apply_flux(Boundary_face& bf) {
  // fetch data
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  double* gh_f = bf.ghost_face(true);
  double* in_f = bf.inside_face(true);
  // initialize to negative of inside flux
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
  // un-invert normal component of momentum
  reflect_normal(gh_f, bf.surface_normal(), nfq, params.n_dim);
}

void Nonpenetration::apply_advection(Boundary_face& bf) {
  Storage_params params = bf.storage_params();
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  double* nrml = bf.surface_normal();
  int nfq = params.n_qpoint()/params.row_size;
  int nd = params.n_dim;
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> veloc(nd);
    Mat<> n(nd);
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      n(i_dim) = nrml[i_dim*nfq + i_qpoint];
      veloc(i_dim) = in_f[i_dim*nfq + i_qpoint];
    }
    veloc -= 2*veloc.dot(n)*n/n.squaredNorm();
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      gh_f[i_dim*nfq + i_qpoint] = veloc(i_dim);
    }
    for (int i_adv = 0; i_adv < params.n_advection(params.row_size); ++i_adv) {
      gh_f[(nd + i_adv)*nfq + i_qpoint] = in_f[(nd + i_adv)*nfq + i_qpoint];
    }
  }
}

No_slip::No_slip(std::shared_ptr<Thermal_bc> thermal, double coercion) : _coercion{coercion}, _thermal{thermal} {}

double Thermal_equilibrium::ghost_heat_flux(Mat<> state, double) {
  double temp = state(last)*.4/state(state.size() - 2)/constants::specific_gas_air;
  double radiative_flux = emissivity*constants::stefan_boltzmann*math::pow(temp, 4);
  double conductive_flux = heat_transfer_coef*(temp - temperature);
  return radiative_flux + conductive_flux;
}

void No_slip::apply_state(Boundary_face& bf) {
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face(false);
  double* in_f = bf.inside_face(false);
  double* sc = bf.state_cache();
  Array<double> presc = bf.prescribed_data();
  int nd = params.n_dim;
  int nfq = params.n_qpoint()/params.row_size; // number of face quadrature points
  // set ghost state
  // momentum
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      gh_f[i_dim*nfq + i_qpoint] = 2*presc(i_dim)[i_qpoint]*in_f[nd*nfq + i_qpoint] - in_f[i_dim*nfq + i_qpoint];
    }
  }
  // density
  for (int i_dof = params.n_dim*nfq; i_dof < (nd + 1)*nfq; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  // energy
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    Mat<> state(params.n_dim + 2); // yes, this should be ignoring turbulence variables
    for (int i_var = 0; i_var < params.n_dim + 2; ++i_var) state(i_var) = in_f[i_var*nfq + i_qpoint];
    gh_f[(params.n_dim + 1)*nfq + i_qpoint] = math::pow(_thermal->ghost_energy(state), 2)/state(last);
  }
  // turbulence variables
  //! \todo __Carter:__ Set the correct wall boundary conditions for the turbulence variables.
  //! Right now, they're set to \f$ k = 1 \f$ and \f$ \omega = 0.01 \f$.
  if (params.n_var == params.n_dim + 4) {
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      double spec_turb_kin_ener_wall = 0.;
      double spec_turb_diss_wall = 4e8;
      gh_f[(params.n_dim + 2)*nfq + i_qpoint] = 2*spec_turb_kin_ener_wall*in_f[params.n_dim*nfq + i_qpoint]
                                                - in_f[(params.n_dim + 2)*nfq + i_qpoint];
      #if 1
      gh_f[(params.n_dim + 3)*nfq + i_qpoint] = 2*std::log(spec_turb_diss_wall)*in_f[params.n_dim*nfq + i_qpoint]
                                                - in_f[(params.n_dim + 3)*nfq + i_qpoint];
      #else
      gh_f[(params.n_dim + 3)*nfq + i_qpoint] = in_f[(params.n_dim + 3)*nfq + i_qpoint];
      #endif
      gh_f[(params.n_dim + 4)*nfq + i_qpoint] = in_f[(params.n_dim + 4)*nfq + i_qpoint];
      gh_f[(params.n_dim + 5)*nfq + i_qpoint] = in_f[(params.n_dim + 5)*nfq + i_qpoint];
      gh_f[(params.n_dim + 6)*nfq + i_qpoint] = in_f[(params.n_dim + 6)*nfq + i_qpoint];
    }
  }
  // prime `state_cache` with average state for use in emissivity BC
  for (int i_dof = 0; i_dof < params.n_var*nfq; ++i_dof) {
    sc[i_dof] = (gh_f[i_dof] + in_f[i_dof])/2;
  }
}

void No_slip::apply_flux(Boundary_face& bf) {
  auto params = bf.storage_params();
  int nfq = params.n_qpoint()/params.row_size;
  double* gh_f = bf.ghost_face(true);
  double* in_f = bf.inside_face(true);
  double* sc = bf.state_cache();
  double* nrml = bf.surface_normal();
  // set momentum and mass flux (pretty straightforward)
  for (int i_dof = 0; i_dof < params.n_dim*nfq; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  for (int i_dof = params.n_dim*nfq; i_dof < (params.n_dim + 1)*nfq; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
  // set energy flux depending on thermal boundary condition
  for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
    int i_dof = i_qpoint + (params.n_dim + 1)*nfq;
    double normal = 0;
    for (int i_dim = 0; i_dim < params.n_dim; ++i_dim) {
      double n = nrml[i_dim*nfq + i_qpoint];
      normal += n*n;
    }
    normal = std::sqrt(normal);
    int flux_sign = 2*bf.inside_face_sign() - 1;
    Mat<> state(params.n_var);
    for (int i_var = 0; i_var < params.n_var; ++i_var) state(i_var) = sc[i_var*nfq + i_qpoint];
    double ghost_heat = _thermal->ghost_heat_flux(state, in_f[i_dof]*flux_sign/normal);
    gh_f[i_dof] = _coercion*(normal*flux_sign*ghost_heat - in_f[i_dof]) + in_f[i_dof];
  }
  // set turbulence variables
  for (int i_dof = (params.n_dim + 2)*nfq; i_dof < params.n_var*nfq; ++i_dof) gh_f[i_dof] = in_f[i_dof];
  for (int i_dof = (params.n_dim + 4)*nfq; i_dof < params.n_var*nfq; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
}

void No_slip::apply_advection(Boundary_face& bf) {
  auto params = bf.storage_params();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint()/params.row_size;
  double* in_f = bf.inside_face(false);
  double* gh_f = bf.ghost_face(false);
  // set velocity to 0
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gh_f[i_dim*nq + i_qpoint] = -in_f[i_dim*nq + i_qpoint];
    }
  }
  // don't change advected scalar
  for (int i_adv = 0; i_adv < params.n_advection(params.row_size); ++i_adv) {
    for (int i_qpoint = 0; i_qpoint < nq; ++i_qpoint) {
      gh_f[(nd + i_adv)*nq + i_qpoint] = in_f[(nd + i_adv)*nq + i_qpoint];
    }
  }
}

void No_slip::set_prescribed(Interpreter& inter, Boundary_face& bf) {
  auto sub = inter.make_sub();
  auto params = bf.storage_params();
  const int nd = params.n_dim;
  const int nq = params.n_qpoint()/params.row_size;
  std::array<std::string, 2> surface_properties {"pos", "normal"};
  Array<double> surface({2, nd, nq});
  surface(0) = Array<double>({nd, nq}, bf.surface_position());
  surface(1) = Array<double>({nd, nq}, bf.surface_normal());
  for (int i = 0; i < 2; ++i) {
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      sub.variables->assign(surface_properties[i] + std::to_string(i_dim), surface(i)(i_dim));
    }
    for (int i_dim = nd; i_dim < 3; ++i_dim) {
      sub.variables->assign(surface_properties[i] + std::to_string(i_dim), 0.);
    }
  }
  auto expr = inter.variables->lookup<std::string>("wall_velocity");
  HEXED_ASSERT(expr, "`wall_velocity` must be specified for `No_slip`");
  sub.exec(expr.value());
  Array<double> data {bf.prescribed_data()};
  HEXED_ASSERT(data.shape()[0] == nd && data.size() == nd*nq,
               "`prescribed_data` is not the right shape to hold the velocity");
  for (int i_dim = 0; i_dim < nd; ++i_dim) {
    sub.variables->assign_array(data(i_dim), "velocity" + std::to_string(i_dim));
  }
}

void Copy::apply_state(Boundary_face& bf) {
  copy_state(bf);
}

void Copy::apply_flux(Boundary_face& bf) {
  copy_state(bf);
}

void Copy::apply_advection(Boundary_face& bf) {
  copy_state(bf);
}

void Outflow::apply_state(Boundary_face& bf) {
  copy_state(bf);
}

void Outflow::apply_flux(Boundary_face& bf) {
  // set to negative of inside flux
  auto params = bf.storage_params();
  double* gh_f = bf.ghost_face(true);
  double* in_f = bf.inside_face(true);
  for (int i_dof = 0; i_dof < params.n_dof()/params.row_size; ++i_dof) gh_f[i_dof] = -in_f[i_dof];
}

Expression_bc::Expression_bc(Interpreter& inter, std::string state_expr, std::string flux_expr)
: _inter{inter}, _exprs{state_expr, flux_expr}
{}

void Expression_bc::_apply(Boundary_face& bf, bool is_flux) {
  auto sub = _inter.make_sub();
  auto params = bf.storage_params();
  const int nd = params.n_dim;
  const int nv = params.n_var;
  const int nq = params.n_qpoint()/params.row_size;
  std::array<std::string, 2> surface_properties {"pos", "normal"};
  std::array<std::string, 2> names {"state", "flux"};
  Array<double> surface({2, nd, nq});
  surface(0) = Array<double>({nd, nq}, bf.surface_position());
  surface(1) = Array<double>({nd, nq}, bf.surface_normal());
  for (int i = 0; i < 2; ++i) {
    for (int i_dim = 0; i_dim < nd; ++i_dim) {
      sub.variables->assign(surface_properties[i] + std::to_string(i_dim), surface(i)(i_dim));
    }
    for (int i_dim = nd; i_dim < 3; ++i_dim) {
      sub.variables->assign(surface_properties[i] + std::to_string(i_dim), 0.);
    }
    Array<double> inside({nv, nq}, bf.inside_face(i));
    for (int i_var = 0; i_var < nv; ++i_var) {
      sub.variables->assign(names[i] + std::to_string(i_var), inside(i_var).copy());
    }
  }
  sub.exec(_exprs[is_flux]);
  Array<double> ghost({nv, nq}, bf.ghost_face(is_flux));
  for (int i_var = 0; i_var < nv; ++i_var) {
    sub.variables->assign_array(ghost(i_var), "ghost_" + names[is_flux] + std::to_string(i_var));
  }
}

}
