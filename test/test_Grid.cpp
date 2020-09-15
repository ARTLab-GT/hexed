#include "catch.hpp"

#include "Grid.hpp"
#include "Equidistant.hpp"

TEST_CASE("Grid")
{
  Equidistant basis (8);
  Grid grid1 (4, 1, 5, 0.1, basis);
  Grid grid2 (4, 2, 5, 0.1, basis);
  Grid grid3 (4, 3, 27, 0.1, basis);
  SECTION("Construction")
  {
    REQUIRE(grid2.basis.node(1) == 1./7.);
    REQUIRE(grid3.n_var == 4);
    REQUIRE(grid1.n_dim == 1);
    REQUIRE(grid2.n_dim == 2);
    REQUIRE(grid3.n_dim == 3);
    REQUIRE(grid1.n_qpoint == 8);
    REQUIRE(grid2.n_qpoint == 64);
    REQUIRE(grid3.n_qpoint == 512);
    REQUIRE(grid1.n_dof == 32);
    REQUIRE(grid2.n_dof == 256);
    REQUIRE(grid3.n_dof == 2048);
    REQUIRE(grid1.n_elem == 5);
    int size = 2048*27;
    REQUIRE(grid3.state_r[0] == 0.);
    REQUIRE(grid3.state_w[0] == 0.);
    REQUIRE(grid3.state_r[size - 1] == 0.);
    REQUIRE(grid3.state_w[size - 1] == 0.);
  }

  for (int i = 0; i < 5; ++i)
  {
    grid1.pos[i] = i - 1;
  }
  {
    int i = 0;
    grid2.pos[i++] = 0; grid2.pos[i++] =  0;
    grid2.pos[i++] = 0; grid2.pos[i++] = -1;
    grid2.pos[i++] = 1; grid2.pos[i++] = -1;
    grid2.pos[i++] = 1; grid2.pos[i++] =  1;
    grid2.pos[i++] = 3; grid2.pos[i++] =  0;
  }
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_qpoint = k + 3*(j + 3*i);
        grid3.pos[i_qpoint + 0] = i;
        grid3.pos[i_qpoint + 1] = j;
        grid3.pos[i_qpoint + 2] = k;
      }
    }
  }

  SECTION("Positioning")
  {
    std::vector<double> pos;
    pos = grid1.get_pos(0);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[1] == 1./7.);
    REQUIRE(pos[7] == 1.);
    pos = grid1.get_pos(1);
    REQUIRE(pos[0] == 1.);

    pos = grid2.get_pos(0);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[1] == 1./7.);
    REQUIRE(pos[8] == 0.);
    REQUIRE(pos[63] == 1.);
    REQUIRE(pos[64] == 0.);
    REQUIRE(pos[65] == 0.);
    REQUIRE(pos[72] == 1./7.);
    REQUIRE(pos[127] == 1.);
    pos = grid2.get_pos(1);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[63] == 1.);
    REQUIRE(pos[64] == -1.);
    REQUIRE(pos[127] == 0.);

    pos = grid3.get_pos(0);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[511] == 1.);
    REQUIRE(pos[512] == 0.);
  }

}