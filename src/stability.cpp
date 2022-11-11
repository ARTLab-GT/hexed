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
};

constexpr int n_side = 10;

int main()
{
  Stability_solver sol(2, 6, 1.);
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
  sol.initialize(hexed::Random_func({0., 0., 1., 2e5}, {0., 0., .1, 2e4}));
  sol.set_art_visc_constant(1.);
  hexed::State_variables state;
  hexed::Component comp(state, 2);
  hexed::Constant_func mean({0., 0., 1., 2e5});
  hexed::Diff_sq diff(state, mean);
  printf("%e\n", sol.bounds_field(diff)[2][1]);
  sol.run_diffusive(1.3e-3);
  printf("%e\n", sol.bounds_field(diff)[2][1]);
  #if HEXED_USE_OTTER
  otter::plot plt;
  sol.visualize_field_otter(plt, comp);
  sol.visualize_edges_otter(plt);
  plt.show();
  #endif
}
