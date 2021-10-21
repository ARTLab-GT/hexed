#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Gauss_legendre.hpp>
#include <Tecplot_file.hpp>

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
  grid.purge_vertices();
  grid.calc_jacobian();
  cartdg::Tecplot_file file {"deformed_test", 2, 4, 0.};
  grid.visualize_interior(file);
}
