#include <nlopt.hpp>
#include <hexed/Solver.hpp>
#include <hexed/kernel_factory.hpp>
#include <hexed/Write_face.hpp>
#include <hexed/Prolong_refined.hpp>

class Stability_solver : public hexed::Solver
{
  public:
  Stability_solver(int n_dim, int row_size, double root_mesh_size)
  : hexed::Solver(n_dim, row_size, root_mesh_size)
  {}
  void run_diffusive(double dt)
  {
    const int nd = params.n_dim;
    const int rs = params.row_size;
    const int n_dof = params.n_dof();
    auto& elems = acc_mesh.elements();
    for (int iter = 0; iter < 100; ++iter) {
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] = state[i_dof];
      }
      update_art_visc(1., true);
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] += dt*state[i_dof];
      }
      (*hexed::kernel_factory<hexed::Write_face>(nd, params.row_size, basis))(elems);
      (*hexed::kernel_factory<hexed::Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
      update_art_visc(1., true, false);
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] += coef[0]*dt*state[i_dof];
      }
      (*hexed::kernel_factory<hexed::Write_face>(nd, params.row_size, basis))(elems);
      (*hexed::kernel_factory<hexed::Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
      update_art_visc(1., true, false);
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof + n_dof] += coef[1]*dt*state[i_dof];
      }
      for (int i_elem = 0; i_elem < elems.size(); ++i_elem) {
        double* state = elems[i_elem].stage(0);
        for (int i_dof = 0; i_dof < n_dof; ++i_dof) state[i_dof] = state[i_dof + n_dof];
      }
      (*hexed::kernel_factory<hexed::Write_face>(nd, params.row_size, basis))(elems);
      (*hexed::kernel_factory<hexed::Prolong_refined>(nd, rs, basis))(acc_mesh.refined_faces());
    }
  }
} sol(2, 6, 1);

hexed::State_variables state;
hexed::Component comp(state, 2);
hexed::Constant_func mean({0., 0., 1., 2e5});
hexed::Diff_sq diff(state, mean);

double objective(const std::vector<double> &x, std::vector<double> &grad, void*)
{
  auto bounds_before = sol.bounds_field(state);
  sol.initialize(hexed::Random_func({0., 0., 1., 2e5}, {0., 0., .1, 2e4}));
  sol.run_diffusive(x[0]);
  auto bounds_after = sol.bounds_field(state);
  double spread_diff =   (bounds_before[2][1] - bounds_before[2][0])
                       - (bounds_after [2][1] - bounds_after [2][0]);
  return spread_diff*spread_diff - x[0]*x[0];
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
  srand(406);
  sol.set_art_visc_constant(1.);
  nlopt::opt opt(nlopt::LN_SBPLX, 1);
  opt.set_min_objective(objective, NULL);
  sol.coef[0] = 6e-4;
  std::vector<double> x {1e-3};
  double min_obj;
  opt.optimize(x, min_obj);
  printf("%e %e\n", x[0], objective(x, x, nullptr));
  #if HEXED_USE_OTTER
  otter::plot plt;
  sol.visualize_field_otter(plt, comp);
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
