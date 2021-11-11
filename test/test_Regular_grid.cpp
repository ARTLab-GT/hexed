#include <catch2/catch.hpp>
#include <iostream>

#include <cartdgConfig.hpp>
#include <Regular_grid.hpp>
#include <Equidistant.hpp>
#include <Gauss_lobatto.hpp>

class Arbitrary_integrand : public cartdg::Domain_func
{
  public:
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state)
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0.};
  }
};

TEST_CASE("Regular_grid")
{
  REQUIRE(cartdg::config::max_row_size >= 8); // too lazy to make this work for other max row size
  const int row_size {8};
  cartdg::Equidistant basis (row_size);
  cartdg::Regular_grid grid1 (4, 1, 5, 0.1, basis);
  cartdg::Regular_grid grid2 (4, 2, 6, 0.1, basis);
  cartdg::Regular_grid grid3 (4, 3, 27, 0.1, basis);
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
    REQUIRE(grid1.time == 0.);
    REQUIRE(grid1.origin.size() == 1);
    REQUIRE(grid2.origin.size() == 2);
    REQUIRE(grid3.origin.size() == 3);
    REQUIRE(grid1.origin[0] == 0.);
    REQUIRE(grid2.origin[0] == 0.);
    REQUIRE(grid2.origin[1] == 0.);
    REQUIRE(grid3.origin[0] == 0.);
    REQUIRE(grid3.origin[1] == 0.);
    REQUIRE(grid3.origin[2] == 0.);
  }

  for (int i = 0; i < 5; ++i)
  {
    grid1.pos[i] = i - 1;
  }
  {
    int i = 0;
    grid2.pos[i++] = 1; grid2.pos[i++] =  1;
    grid2.pos[i++] = 0; grid2.pos[i++] =  0;
    grid2.pos[i++] = 0; grid2.pos[i++] = -1;
    grid2.pos[i++] = 1; grid2.pos[i++] = -1;
    grid2.pos[i++] = 3; grid2.pos[i++] =  0;
    grid2.pos[i++] = -1; grid2.pos[i++] = -1;
  }
  grid2.origin[1] = 10.;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_elem = 3*(k + 3*(j + 3*i));
        grid3.pos[i_elem + 0] = i;
        grid3.pos[i_elem + 1] = j;
        grid3.pos[i_elem + 2] = k;
      }
    }
  }

  SECTION("Positioning")
  {
    std::vector<double> pos;

    pos = grid1.get_pos(0);
    REQUIRE(pos[0] == -0.1);
    REQUIRE(pos[1] == Approx(-0.1 + 0.1/7.));
    REQUIRE(pos[7] == 0.);
    pos = grid1.get_pos(1);
    REQUIRE(pos[0] == 0.);

    pos = grid2.get_pos(0);
    REQUIRE(pos[0] == 0.1);
    pos = grid2.get_pos(1);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[1] == 0.);
    REQUIRE(pos[8] == Approx(0.1/7.));
    REQUIRE(pos[63] == 0.1);
    REQUIRE(pos[64] == 10.);
    REQUIRE(pos[65] == Approx(10. + 0.1/7.));
    REQUIRE(pos[72] == Approx(10.));
    REQUIRE(pos[127] == 10.1);
    pos = grid2.get_pos(2);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[63] == 0.1);
    REQUIRE(pos[64] == 10. - 0.1);
    REQUIRE(pos[127] == 10.);

    pos = grid3.get_pos(0);
    REQUIRE(pos[0] == 0.);
    REQUIRE(pos[511] == 0.1);
    REQUIRE(pos[512] == 0.);
    pos = grid3.get_pos(1);
    REQUIRE(pos[0  ] == 0.);
    REQUIRE(pos[512] == 0.);
    REQUIRE(pos[1024] == 0.1);
    pos = grid3.get_pos(3);
    REQUIRE(pos[0  ] == 0.);
    REQUIRE(pos[512] == 0.1);
    REQUIRE(pos[1024] == 0.);
    pos = grid3.get_pos(9);
    REQUIRE(pos[0  ] == 0.1);
    REQUIRE(pos[512] == 0.);
    REQUIRE(pos[1024] == 0.);
    pos = grid3.get_pos(26);
    REQUIRE(pos[0  ] == 0.2);
    REQUIRE(pos[512] == 0.2);
    REQUIRE(pos[1024] == 0.2);
  }

  SECTION("State integration")
  {
    cartdg::Gauss_lobatto basis (5);
    cartdg::Regular_grid grid1 (2, 1, 5, 0.1, basis);
    cartdg::Regular_grid grid2 (2, 2, 5, 0.1, basis);
    cartdg::Regular_grid grid3 (2, 3, 5, 0.1, basis);
    for (cartdg::Regular_grid* grid : {&grid1, &grid2, &grid3})
    {
      for (int i_elem = 0; i_elem < grid->n_elem; ++i_elem)
      {
        double* stage = grid->element(i_elem).stage(0);
        for (int i_dof = 0; i_dof < grid->n_dof; ++i_dof)
        {
          stage[i_dof] = 1.;
        }
      }
    }
    REQUIRE(grid1.integral()[0] == Approx(0.5  ).margin(1e-12));
    REQUIRE(grid1.integral()[1] == Approx(0.5  ).margin(1e-12));
    REQUIRE(grid2.integral()[0] == Approx(0.05 ).margin(1e-12));
    REQUIRE(grid2.integral()[1] == Approx(0.05 ).margin(1e-12));
    REQUIRE(grid3.integral()[0] == Approx(0.005).margin(1e-12));
    REQUIRE(grid3.integral()[1] == Approx(0.005).margin(1e-12));

    cartdg::Regular_grid square (1, 2, 2, 1., basis);
    for (int i_elem = 0; i_elem < square.n_elem; ++i_elem)
    {
      double* stage = square.element(i_elem).stage(0);
      for (int i = 0; i < basis.row_size; ++i)
      {
        for (int j = 0; j < basis.row_size; ++j)
        {
          double pos0 = basis.node(i) + i_elem; double pos1 = basis.node(j);
          stage[i*basis.row_size + j] = pos0*pos0*pos1*pos1*pos1;
          square.pos[i_elem*2] = i_elem;
        }
      }
    }
    square.time = 0.61;
    REQUIRE(square.integral()[0] == Approx(2./3.));
    Arbitrary_integrand integrand;
    REQUIRE(square.integral(integrand)[0] == Approx(2*0.61));
  }

  SECTION("Automatic graph creation")
  {
    SECTION("non-periodic")
    {
      grid2.auto_connect();
      REQUIRE(grid2.n_con(0) == 2);
      REQUIRE(grid2.n_con(1) == 1);
      {
        auto con = grid2.connection(0, 0);
        REQUIRE(con[0] == &grid2.element(2));
        REQUIRE(con[1] == &grid2.element(3));
      }
      {
        auto con = grid2.connection(0, 1);
        REQUIRE(con[0] == &grid2.element(5));
        REQUIRE(con[1] == &grid2.element(2));
      }
      {
        auto con = grid2.connection(1, 0);
        REQUIRE(con[0] == &grid2.element(2));
        REQUIRE(con[1] == &grid2.element(1));
      }
    }

    SECTION("periodic")
    {
      std::vector<int> periods {0, 3};
      grid2.auto_connect(periods);
      REQUIRE(grid2.n_con(0) == 2);
      REQUIRE(grid2.n_con(1) == 2);

      std::vector<int> periods3d {3, 3, 3};
      grid3.auto_connect(periods3d);
      REQUIRE(grid3.n_con(0) == 27);
      REQUIRE(grid3.n_con(1) == 27);
      REQUIRE(grid3.n_con(2) == 27);
    }
  }

  SECTION("Runge Kutta time integration")
  {
    for (int i_elem = 0; i_elem < grid1.n_elem; ++i_elem)
    {
      for (int i_dof = 0; i_dof < grid1.n_dof; ++i_dof)
      {
        grid1.element(i_elem).stage(grid1.i_stage_read())[i_dof] = 7.;
      }
    }

    do
    {
      for (int i_elem = 0; i_elem < grid1.n_elem; ++i_elem)
      {
        for (int i_dof = 0; i_dof < grid1.n_dof; ++i_dof)
        {
          double* read = grid1.element(i_elem).stage(grid1.i_stage_read());
          double* write = grid1.element(i_elem).stage(grid1.i_stage_write());
          write[i_dof] = read[i_dof] + 0.371;
        }
      }
    }
    while (!grid1.execute_runge_kutta_stage());

    for (int i_elem = 0; i_elem < grid1.n_elem; ++i_elem)
    {
      for (int i_dof = 0; i_dof < grid1.n_dof; ++i_dof)
      {
        REQUIRE(grid1.element(i_elem).stage(grid1.i_stage_read())[i_dof] == Approx(7.371));
      }
    }

    do
    {
      for (int i_elem = 0; i_elem < grid1.n_elem; ++i_elem)
      {
        for (int i_dof = 0; i_dof < grid1.n_dof; ++i_dof)
        {
          double* read = grid1.element(i_elem).stage(grid1.i_stage_read());
          double* write = grid1.element(i_elem).stage(grid1.i_stage_write());
          write[i_dof] = read[i_dof] + 0.001;
        }
      }
    }
    while (!grid1.execute_runge_kutta_stage());

    for (int i_elem = 0; i_elem < grid1.n_elem; ++i_elem)
    {
      for (int i_dof = 0; i_dof < grid1.n_dof; ++i_dof)
      {
        REQUIRE(grid1.element(i_elem).stage(grid1.i_stage_read())[i_dof] == Approx(7.372));
      }
    }
  }

  SECTION("Add element")
  {
    cartdg::Equidistant basis (8);
    cartdg::Regular_grid grid (4, 1, 0, 0.1, basis);
    std::vector<int> position {1};
    REQUIRE(grid.n_elem == 0);
    REQUIRE(grid.add_element(position) == 0);
    REQUIRE(grid.n_elem == 1);
    grid.element(0).stage(0)[0] = 1.; // modify element to verify existence
    REQUIRE(grid.element(0).stage(0)[0] == 1.);
    position[0] += 1;
    REQUIRE(grid.add_element(position) == 1);
    REQUIRE(grid.n_elem == 2);
    REQUIRE(grid.pos[0] == 1);
    REQUIRE(grid.pos[1] == 2);
    grid.element(1).stage(0)[0] = 1.;
    REQUIRE(grid.element(1).stage(0)[0] == 1.);
    grid.auto_connect();
    REQUIRE(grid.n_con(0) == 1);
  }
}
