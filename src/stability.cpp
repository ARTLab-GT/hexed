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

  void init()
  {
    srand(406);
    std::vector<double> init(params.n_dim, 0.);
    init.push_back(1.);
    init.push_back(2e-5);
    std::vector<double> var;
    for (double elem : init) var.push_back(.1*elem);
    initialize(hexed::Random_func(init, var));
  }

  double amplification(double dt)
  {
    init();
    const int nd = params.n_dim;
    auto bounds_before = bounds_field(state);
    run_diffusive(dt);
    auto bounds_after = bounds_field(state);
    return   (bounds_after [nd][1] - bounds_after [nd][0])
           - (bounds_before[nd][1] - bounds_before[nd][0]);
  }

  double time_step()
  {
    try {
      double result = hexed::custom_math::bisection([this](double dt){return amplification(dt);}, {1e-6, 1e-1}, 1e-7);
      return (result < .1) ? result : 0;
    } catch (const std::runtime_error& e) {
      std::cerr << "`time_step` returned 0 due to following exception:\n  " << e.what() << "\n";
      return 0;
    }
  }
};

const int n_dim = 1;
Stability_solver sol(n_dim, 6, .13);

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

int main()
{
  int n_side = 5;
  int n_elem = hexed::custom_math::pow(n_side, n_dim);
  std::vector<int> sn(n_elem);
  std::vector<int> strides(n_dim);
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
    strides[i_dim] = hexed::custom_math::pow(n_side, n_dim - 1 - i_dim);
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    std::vector<int> pos(n_dim);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      pos[i_dim] = (i_elem/strides[i_dim])%n_side;
    }
    sn[i_elem] = sol.mesh().add_element(0, 0, pos);
  }
  for (int i_elem = 0; i_elem < n_elem; ++i_elem) {
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      int row = (sn[i_elem]/strides[i_dim])%n_side;
      int new_row = (row + 1)%n_side;
      int new_sn = sn[i_elem] + (new_row - row)*strides[i_dim];
      sol.mesh().connect_cartesian(0, {sn[i_elem], new_sn}, {i_dim});
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
  sol.init();
  for (int i = 0; i < 1000; ++i) sol.update(.9);
  #if HEXED_USE_OTTER
  otter::plot plt;
  sol.visualize_field_otter(plt, comp);
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
