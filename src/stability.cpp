#include <nlopt.hpp>
#include <hexed/Solver.hpp>
#include <hexed/kernel_factory.hpp>
#include <hexed/Write_face.hpp>
#include <hexed/Prolong_refined.hpp>
#include <hexed/Mcs_cartesian.hpp>
#include <iostream>

hexed::State_variables state;
hexed::Component comp(state, 2);
hexed::Constant_func mean({0., 0., 1., 2e5});
hexed::Diff_sq diff(state, mean);

class Stability_solver : public hexed::Solver
{
  public:
  bool square = false;
  double coef [2] {};
  Stability_solver(int n_dim, int row_size, double root_mesh_size)
  : hexed::Solver(n_dim, row_size, root_mesh_size)
  {}
  void run_diffusive(double cfl)
  {
    const int nd = params.n_dim;
    const int rs = params.row_size;
    const int n_dof = params.n_dof();
    auto& elems = acc_mesh.elements();
    double mcs = (*hexed::kernel_factory<hexed::Mcs_cartesian>(nd, rs, hexed::char_speed::Art_visc(), 2, 1))(acc_mesh.cartesian().elements());
    double dt = cfl/mcs;
    for (int iter = 0; iter < 100; ++iter) {
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] = state[i_dof];
      }
      update_laplacian(1., true);
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] += dt*state[i_dof];
      }
      update_laplacian(1., true, false);
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) {
          state[i_dof] = coef[0]*dt*(square ? cfl/coef[1] : 1.)/mcs*state[i_dof] + state[i_dof + n_dof];
        }
      }
      (*hexed::kernel_factory<hexed::Write_face>(nd, params.row_size, basis))(elems);
      (*hexed::kernel_factory<hexed::Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
    }
  }

  double amplification(double dt)
  {
    srand(406);
    initialize(hexed::Random_func({0., 0., 1., 2e5}, {0., 0., .1, 2e4}));
    auto bounds_before = bounds_field(state);
    run_diffusive(dt);
    auto bounds_after = bounds_field(state);
    return   (bounds_after [2][1] - bounds_after [2][0])
           - (bounds_before[2][1] - bounds_before[2][0]);
  }

  double time_step()
  {
    double result = hexed::custom_math::bisection([this](double dt){return amplification(dt);}, {1e-7, 1}, 1e-7);
    return (result < .1) ? result : 0;
  }
} sol(2, 6, .13);

int i = 0;

double objective(const std::vector<double> &x, std::vector<double> &grad, void*)
{
  for (int j = 0; j < int(x.size()); ++j) {
    printf("%e ", x[j]);
    sol.coef[i + j] = x[j];
  }
  std::cout << std::flush;
  double dt = sol.time_step();
  printf("| %e\n", dt);
  return -dt;
}

constexpr int n_side = 10;

int main()
{
  int sn [n_side][n_side];
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      sn[i][j] = sol.mesh().add_element(0, 0, {i, j});
    }
  }
  for (int i = 0; i < n_side; ++i) {
    for (int j = 0; j < n_side; ++j) {
      sol.mesh().connect_cartesian(0, {sn[i][j], sn[(i + 1)%n_side][j]}, {0});
      sol.mesh().connect_cartesian(0, {sn[i][j], sn[i][(j + 1)%n_side]}, {1});
    }
  }
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.set_local_tss();
  sol.set_art_visc_constant(3.4);
  printf("unmodified time step: %e\n", sol.time_step());
  nlopt::opt opt(nlopt::LN_SBPLX, 1);
  opt.set_min_objective(objective, nullptr);
  opt.set_xtol_rel(1e-2);
  std::vector<double> x {1e-5};
  double min_obj;
  opt.optimize(x, min_obj);
  double cfl = -objective(x, x, nullptr);
  sol.coef[0] = x[0];
  sol.coef[1] = cfl*.95;
  printf("cancellation: {%e}\nCFL: %e\n", sol.coef[0], sol.coef[1]);
  int n = 20;
  sol.square = true;
  for (int i = 0; i <= n; ++i) {
    double stab_rat = double(i)/n;
    printf("%e: %e\n", stab_rat, sol.amplification(stab_rat*sol.coef[1]));
  }
  srand(406);
  sol.initialize(hexed::Random_func({0., 0., 1., 2e5}, {0., 0., .1, 2e4}));
  for (int i = 0; i < 1000; ++i) sol.update(.9);
  #if HEXED_USE_OTTER
  otter::plot plt;
  sol.visualize_field_otter(plt, comp);
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
