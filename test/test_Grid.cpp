#include <catch2/catch.hpp>
#include <iostream>

#include <cartdgConfig.hpp>
#include <Grid.hpp>
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

TEST_CASE("Grid")
{
  cartdg::Equidistant basis (CARTDG_MAX_BASIS_ROW_SIZE);
  cartdg::Grid grid1 (4, 1, 5, 0.1, basis);
  cartdg::Grid grid2 (4, 2, 6, 0.1, basis);
  cartdg::Grid grid3 (4, 3, 27, 0.1, basis);
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

    int size = 2048*27;
    REQUIRE(grid3.state_r()[0] == 0.);
    REQUIRE(grid3.state_w()[0] == 0.);
    REQUIRE(grid3.state_r()[size - 1] == 0.);
    REQUIRE(grid3.state_w()[size - 1] == 0.);
    size = 27*std::pow(CARTDG_MAX_BASIS_ROW_SIZE, 3);
    REQUIRE((int)grid3.derivs.size() == size);
    REQUIRE(grid3.derivs[0] == 0.);
    REQUIRE(grid3.derivs[size - 1] == 0.);
    REQUIRE((int)grid3.visc.size() == 27*8);
    REQUIRE(grid3.visc[0] == 0.);
    REQUIRE(grid3.visc[27*8 - 1] == 0.);
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

  SECTION("(trivial) Jacobian determinant")
  {
    REQUIRE(grid1.jacobian_det(0, 0) == 1.);
    REQUIRE(grid2.jacobian_det(4, 2) == 1.);
    REQUIRE(grid3.jacobian_det(16, 8) == 1.);
    REQUIRE(grid3.jacobian_det(0, 0) == 1.);
  }

  SECTION("Visualization")
  {
    grid1.visualize("unit_test_1d");
    grid2.visualize("unit_test_2d");
    grid3.visualize("unit_test_3d");
  }

  SECTION("State integration")
  {
    cartdg::Gauss_lobatto basis (5);
    cartdg::Grid grid1 (2, 1, 5, 0.1, basis);
    cartdg::Grid grid2 (2, 2, 5, 0.1, basis);
    cartdg::Grid grid3 (2, 3, 5, 0.1, basis);
    for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i) grid1.state_r()[i] = 1.;
    for (int i = 0; i < grid2.n_dof*grid1.n_elem; ++i) grid2.state_r()[i] = 1.;
    for (int i = 0; i < grid3.n_dof*grid1.n_elem; ++i) grid3.state_r()[i] = 1.;
    REQUIRE(grid1.integral()[0] == Approx(0.5  ).margin(1e-12));
    REQUIRE(grid1.integral()[1] == Approx(0.5  ).margin(1e-12));
    REQUIRE(grid2.integral()[0] == Approx(0.05 ).margin(1e-12));
    REQUIRE(grid2.integral()[1] == Approx(0.05 ).margin(1e-12));
    REQUIRE(grid3.integral()[0] == Approx(0.005).margin(1e-12));
    REQUIRE(grid3.integral()[1] == Approx(0.005).margin(1e-12));

    cartdg::Grid square (1, 2, 2, 1., basis);
    for (int i_elem = 0; i_elem < square.n_elem; ++i_elem)
    {
      for (int i = 0; i < basis.row_size; ++i)
      {
        for (int j = 0; j < basis.row_size; ++j)
        {
          double pos0 = basis.node(i) + i_elem; double pos1 = basis.node(j);
          square.state_r()[(i_elem*basis.row_size + i)*basis.row_size + j] = pos0*pos0*pos1*pos1*pos1;
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
    grid2.auto_connect();
    std::vector<int> n_neighb_con = grid2.n_neighb_con();
    REQUIRE(n_neighb_con[0] == 2);
    REQUIRE(n_neighb_con[1] == 1);
    REQUIRE(grid2.neighbor_connections_r()[0][0] == grid2.state_r() + 2*grid2.n_dof);
    REQUIRE(grid2.neighbor_connections_r()[0][1] == grid2.state_r() + 3*grid2.n_dof);
    REQUIRE(grid2.neighbor_connections_r()[0][2] == grid2.state_r() + 5*grid2.n_dof);
    REQUIRE(grid2.neighbor_connections_r()[0][3] == grid2.state_r() + 2*grid2.n_dof);
    REQUIRE(grid2.deriv_neighbor_connections()[0][0] == grid2.derivs.data() + 2*grid2.n_qpoint);
    REQUIRE(grid2.deriv_neighbor_connections()[0][1] == grid2.derivs.data() + 3*grid2.n_qpoint);
    REQUIRE(grid2.deriv_neighbor_connections()[0][2] == grid2.derivs.data() + 5*grid2.n_qpoint);
    REQUIRE(grid2.deriv_neighbor_connections()[0][3] == grid2.derivs.data() + 2*grid2.n_qpoint);
    REQUIRE(grid2.visc_neighbor_connections()[0][0] == grid2.visc.data() + 2*4);
    REQUIRE(grid2.visc_neighbor_connections()[0][1] == grid2.visc.data() + 3*4);
    REQUIRE(grid2.visc_neighbor_connections()[0][2] == grid2.visc.data() + 5*4);
    REQUIRE(grid2.visc_neighbor_connections()[0][3] == grid2.visc.data() + 2*4);
    grid2.clear_neighbors();
    std::vector<int> periods {0, 3};
    grid2.auto_connect(periods);
    n_neighb_con = grid2.n_neighb_con();
    REQUIRE(n_neighb_con[0] == 2);
    REQUIRE(n_neighb_con[1] == 2);

    std::vector<int> periods3d {3, 3, 3};
    grid3.auto_connect(periods3d);
    n_neighb_con = grid3.n_neighb_con();
    REQUIRE(n_neighb_con[0] == 27);
    REQUIRE(n_neighb_con[1] == 27);
    REQUIRE(n_neighb_con[2] == 27);
  }

  SECTION("Runge Kutta time integration")
  {
    for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i)
    {
      grid1.state_r()[i] = 7.;
    }

    do
    {
      for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i)
      {
        grid1.state_w()[i] = grid1.state_r()[i] + 0.371;
      }
    }
    while (!grid1.execute_runge_kutta_stage());

    for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i)
    {
      REQUIRE(grid1.state_r()[i] == Approx(7.371));
    }

    do
    {
      for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i)
      {
        grid1.state_w()[i] = grid1.state_r()[i] + 0.001;
      }
    }
    while (!grid1.execute_runge_kutta_stage());

    for (int i = 0; i < grid1.n_dof*grid1.n_elem; ++i)
    {
      REQUIRE(grid1.state_r()[i] == Approx(7.372));
    }
  }

  SECTION("Add element")
  {
    cartdg::Equidistant basis (8);
    cartdg::Grid grid (4, 1, 0, 0.1, basis);
    std::vector<int> position {1};
    REQUIRE(grid.state_storage[0].size() == 0);
    REQUIRE(grid.derivs.size() == 0);
    REQUIRE(grid.visc.size() == 0);
    REQUIRE(grid.add_element(position) == 0);
    REQUIRE(grid.state_storage[0].size() == 4*8);
    REQUIRE(grid.derivs.size() == 8);
    REQUIRE(grid.visc.size() == 2);
    position[0] += 1;
    REQUIRE(grid.add_element(position) == 1);
    REQUIRE(grid.n_elem == 2);
    REQUIRE(grid.pos[0] == 1);
    REQUIRE(grid.pos[1] == 2);
    grid.auto_connect();
    REQUIRE(grid.n_neighb_con()[0] == 1);

    grid.state_w()[0] = 1.;
    REQUIRE(grid.state_r()[0] == 0.);
    grid.state_w()[2*grid.n_dof - 1] = 1.;
    REQUIRE(grid.state_r()[2*grid.n_dof - 1] == 0.);
  }

}
