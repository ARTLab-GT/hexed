#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Gauss_legendre.hpp>
#include <Tecplot_file.hpp>

class Copy_bc : public cartdg::Ghost_boundary_condition
{
  public:
  Copy_bc(cartdg::Grid& grid, int id, bool ip)
  : cartdg::Ghost_boundary_condition{grid, id, ip}
  {}
  virtual void calc_ghost_state()
  {
    ghost_state() = domain_state();
  }
};

TEST_CASE("Deformed elements")
{
  const int row_size {CARTDG_MAX_BASIS_ROW_SIZE};
  const int n_qpoint {row_size*row_size};
  cartdg::Gauss_legendre basis {row_size};
  cartdg::Deformed_grid grid {4, 2, 0, 0.2, basis};
  double veloc [2] {1.2, -1.3};
  double pres {101000};

  for (int i = 0, i_elem = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      grid.add_element({i, j});
      if (i > 0) grid.connect({i_elem - 3, i_elem}, {0, 0}, {1, 0});
      if (j > 0) grid.connect({i_elem - 1, i_elem}, {1, 1}, {1, 0});
      auto pos = grid.get_pos(i_elem);
      double* stage {grid.element(i_elem).stage(0)};
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        double mass = 1. + 0.03*pos[i_qpoint] + 0.02*pos[i_qpoint + n_qpoint];
        stage[0*n_qpoint + i_qpoint] = mass*veloc[0];
        stage[1*n_qpoint + i_qpoint] = mass*veloc[1];
        stage[2*n_qpoint + i_qpoint] = mass;
        stage[3*n_qpoint + i_qpoint] = pres/0.4 + 0.5*mass*(veloc[0]*veloc[0] + veloc[1]*veloc[1]);
      }
      ++i_elem;
    }
  }

  std::vector<Copy_bc> bcs;
  for (int i_dim : {0, 1})
  {
    for (bool is_positive : {0, 1})
    {
      bcs.emplace_back(grid, i_dim, is_positive);
      for (int i_row = 0; i_row < 3; ++i_row)
      {
        int stride {i_dim ? 1 : 3};
        bcs.back().add_element(i_row*3/stride + is_positive*2*stride);
      }
    }
  }
  for (Copy_bc& bc : bcs) grid.ghost_bound_conds.push_back(&bc);
  grid.purge_vertices();
  grid.calc_jacobian();
  cartdg::Kernel_settings settings;
  double dt = 0.07;
  settings.d_t_by_d_pos = 0.07/0.2;

  SECTION("plain")
  {
    grid.execute_write_face(settings);
    grid.execute_neighbor(settings);
    grid.execute_local(settings);
    for (int i_elem = 0; i_elem < 9; ++i_elem)
    {
      double* r {grid.element(i_elem).stage(settings.i_read)};
      double* w {grid.element(i_elem).stage(settings.i_write)};
      for (int i_dof = 0; i_dof < 4*n_qpoint; ++i_dof)
      {
        r[i_dof] = (w[i_dof] - r[i_dof])/dt;
      }
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        REQUIRE(r[2*n_qpoint + i_qpoint] == Approx(-0.03*veloc[0] - 0.02*veloc[1]));
      }
    }
    cartdg::Tecplot_file file {"deformed_test", 2, 4, 0.};
    grid.visualize_interior(file);
  }
}
