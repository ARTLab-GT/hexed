#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Regular_grid.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("Deformed grid class")
{
  const int row_size = std::min<int>(3, CARTDG_MAX_BASIS_ROW_SIZE);
  cartdg::Gauss_lobatto basis (row_size);
  std::vector<cartdg::Deformed_grid> grids;
  grids.emplace_back(1, 2, 0, 0.2, basis);
  grids.emplace_back(1, 3, 0, 0.2, basis);
  cartdg::Deformed_grid& grid2 = grids[0];
  cartdg::Deformed_grid& grid3 = grids[1];

  SECTION("construction")
  {
    for (cartdg::Deformed_grid& grid : grids)
    {
      REQUIRE(grid.i_elem_wall.empty());
      REQUIRE(grid.i_dim_wall.empty());
      REQUIRE(grid.is_positive_wall.empty());
    }
    REQUIRE(grid2.n_vertices == 4);
    REQUIRE(grid3.n_vertices == 8);
    REQUIRE_THROWS(cartdg::Deformed_grid (1, 2, 1, 0.2, basis));
  }

  SECTION("Adding deformed elements")
  {
    REQUIRE(grid2.n_elem == 0);
    grid2.origin[0] = 0.1;
    grid2.add_element(std::vector<int>{-1, 2});
    grid2.add_element(std::vector<int>{ 1, 2});
    REQUIRE(grid2.n_elem == 2);
    REQUIRE(grid2.deformed_element(0).vertex(0).pos[0] == Approx(-0.1));
    REQUIRE(grid2.deformed_element(0).vertex(1).pos[1] == Approx( 0.6));
    REQUIRE(grid2.deformed_element(1).vertex(0).pos[0] == Approx( 0.3));
    // check that deformed_element(int) and element(int) point to the same thing
    grid2.deformed_element(1).stage(0)[0] = -0.31;
    REQUIRE(grid2.element(1).stage(0)[0] == -0.31);

    REQUIRE(grid3.n_elem == 0);
    grid3.add_element(std::vector<int>{2, -3, 1});
    REQUIRE(grid3.n_elem == 1);
    REQUIRE(grid3.deformed_element(0).vertex(0).pos[0] == Approx( 0.4));

    SECTION("vertex interpolation")
    {
      /*
      Approximate vertex layout:
         3
                 2


      1
                   0
      */
      grid2.deformed_element(0).vertex(0).pos = { 0.1, 2.0, 0.0};
      grid2.deformed_element(0).vertex(1).pos = {-1.1, 2.1, 0.0};
      grid2.deformed_element(0).vertex(2).pos = { 0.0, 2.8, 0.0};
      grid2.deformed_element(0).vertex(3).pos = {-1.0, 2.9, 0.0};
      std::vector<double> pos2 = grid2.get_pos(0);
      REQUIRE(pos2[ 0] == Approx( 0.1));
      REQUIRE(pos2[ 1] == Approx(-0.5));
      REQUIRE(pos2[ 2] == Approx(-1.1));
      REQUIRE(pos2[ 4] == Approx(-0.5));
      REQUIRE(pos2[ 9] == Approx( 2. ));
      REQUIRE(pos2[12] == Approx( 2.4));
      REQUIRE(pos2[13] == Approx( 2.45));
      REQUIRE(pos2[17] == Approx( 2.9));

      grid3.add_element({0, 0, 0});
      grid3.deformed_element(1).vertex(7).pos = {0.2*0.8, 0.2*0.8, 0.2*0.8};
      std::vector<double> pos3 = grid3.get_pos(1);
      for (int i_dim = 0; i_dim < 3; ++i_dim) REQUIRE(pos3[27*i_dim + 13] == Approx(0.2*0.475));
    }

    SECTION("warped face interpolation")
    {
      cartdg::Deformed_element& elem = grid2.deformed_element(0);
      elem.vertex(0).pos = {0.0, 0.0, 0.0};
      elem.vertex(1).pos = {0.0, 1.0, 0.0};
      elem.vertex(2).pos = {1.0, 0.0, 0.0};
      elem.vertex(3).pos = {0.6, 1.0, 0.0};
      elem.node_adjustments()[2*3 + 1] =  0.2;
      elem.node_adjustments()[3*3 + 1] = -0.1;
      std::vector<double> pos2 = grid2.get_pos(0);
      REQUIRE(pos2[ 0] == Approx(0.0));
      REQUIRE(pos2[ 7] == Approx(0.8));
      REQUIRE(pos2[ 8] == Approx(0.6));
      REQUIRE(pos2[16] == Approx(0.5));

      REQUIRE(pos2[ 3] == Approx(0.5 - 0.2*0.2));
      REQUIRE(pos2[ 4] == Approx(0.4 - 0.2*(0.2 - 0.1)/2));
      REQUIRE(pos2[ 5] == Approx(0.3 + 0.2*0.1));
      REQUIRE(pos2[12] == Approx(0.0 + 0.2));
      REQUIRE(pos2[13] == Approx(0.5 + (0.2 - 0.1)/2));
      REQUIRE(pos2[14] == Approx(1.0 - 0.1));

      cartdg::Deformed_element& elem1 = grid2.deformed_element(1);
      elem1.vertex(0).pos = {0.0, 0.0, 0.0};
      elem1.vertex(1).pos = {0.0, 1.0, 0.0};
      elem1.vertex(2).pos = {1.0, 0.0, 0.0};
      elem1.vertex(3).pos = {1.0, 1.0, 0.0};
      elem1.node_adjustments()[1] =  0.1;
      std::vector<double> pos21 = grid2.get_pos(1);
      REQUIRE(pos21[ 0] == Approx(0.0));
      REQUIRE(pos21[ 6] == Approx(1.0));
      REQUIRE(pos21[ 4] == Approx(0.55));

      grid3.add_element({0, 0, 0});
      grid3.deformed_element(1).node_adjustments()[4] = 0.01;
      std::vector<double> pos3 = grid3.get_pos(1);
      REQUIRE(pos3[13] == 0.101);
      REQUIRE(pos3[13 + 27] == .1);
      REQUIRE(pos3[13 + 2*27] == .1);
    }
  }

  SECTION("Element connection")
  {
    SECTION("Same direction")
    {
      grid3.add_element({0, 0, 0});
      grid3.add_element({1, 0, 0});
      grid3.connect({0, 1}, {0, 0}, {1, 0});
      REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(1).vertex(0));
      REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(1));
      REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(2));
      REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(3));

      grid3.add_element({0, 1, 0});
      grid3.connect({2, 0}, {1, 1}, {0, 1});
      REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(2).vertex(0));
      REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(2).vertex(1));
      REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(2).vertex(4));
      REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(2).vertex(5));

      grid3.add_element({0, 0, -1});
      grid3.connect({0, 3}, {2, 2}, {0, 1});
      REQUIRE(&grid3.deformed_element(0).vertex(0) == &grid3.deformed_element(3).vertex(1));
      REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(3).vertex(3));
      REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(3).vertex(5));
      REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(3).vertex(7));

      REQUIRE(grid3.connection(1).element[0] == &grid3.deformed_element(2));
      REQUIRE(grid3.connection(1).element[1] == &grid3.deformed_element(0));
    }
    
    SECTION("Different direction")
    {
      SECTION("0+ 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, -1, 0});
        grid3.connect({0, 1}, {0, 1}, {1, 1});
        REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(1).vertex(2));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(6));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(7));
      }

      SECTION("0- 1-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({-1, 1, 0});
        grid3.connect({0, 1}, {0, 1}, {0, 0});
        REQUIRE(&grid3.deformed_element(0).vertex(0) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(1) == &grid3.deformed_element(1).vertex(1));
        REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(1).vertex(4));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(5));
      }

      SECTION("2+ 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({0, -1, 1});
        grid3.connect({0, 1}, {2, 1}, {1, 1});
        REQUIRE(&grid3.deformed_element(0).vertex(1) == &grid3.deformed_element(1).vertex(2));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(6));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(7));
      }

      SECTION("0+ 2+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, -1});
        grid3.connect({0, 1}, {0, 2}, {1, 1});
        REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(1).vertex(1));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(5));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(7));
      }

      SECTION("2+ 0+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, -1});
        grid3.connect({1, 0}, {2, 0}, {1, 1});
        REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(1).vertex(1));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(5));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(7));
      }

      SECTION("0+ 1-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({0, 1}, {0, 1}, {1, 0});
        REQUIRE(&grid3.deformed_element(0).vertex(4) == &grid3.deformed_element(1).vertex(4));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(5));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(1));

        REQUIRE(grid3.connection(0).i_dim[0] == 0);
        REQUIRE(grid3.connection(0).i_dim[1] == 1);
        REQUIRE(grid3.connection(0).is_positive[0] == true);
        REQUIRE(grid3.connection(0).is_positive[1] == false);
      }

      SECTION("0- 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({1, 0}, {0, 1}, {0, 1});
        REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(1).vertex(2));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(1));
      }

      SECTION("1+ 0-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({0, 1}, {1, 0}, {1, 0});
        REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(1).vertex(2));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(1));
      }

      SECTION("1+ 2-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({0, 1, 1});
        grid3.connect({0, 1}, {1, 2}, {1, 0});
        REQUIRE(&grid3.deformed_element(0).vertex(2) == &grid3.deformed_element(1).vertex(2));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(6) == &grid3.deformed_element(1).vertex(6));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(4));
      }

      SECTION("2+ 0-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, 1});
        grid3.connect({0, 1}, {2, 0}, {1, 0});
        REQUIRE(&grid3.deformed_element(0).vertex(1) == &grid3.deformed_element(1).vertex(1));
        REQUIRE(&grid3.deformed_element(0).vertex(3) == &grid3.deformed_element(1).vertex(3));
        REQUIRE(&grid3.deformed_element(0).vertex(5) == &grid3.deformed_element(1).vertex(0));
        REQUIRE(&grid3.deformed_element(0).vertex(7) == &grid3.deformed_element(1).vertex(2));
        grid3.visualize("test");
      }

      SECTION("2D")
      {
        grid2.add_element({0, 0});
        grid2.add_element({1, 0});
        grid2.connect({0, 1}, {0, 0}, {1, 0});
        REQUIRE(&grid2.deformed_element(0).vertex(2) == &grid2.deformed_element(1).vertex(0));
        REQUIRE(&grid2.deformed_element(0).vertex(3) == &grid2.deformed_element(1).vertex(1));

        grid2.add_element({-1, -1});
        grid2.connect({0, 2}, {0, 1}, {0, 1});
        REQUIRE(&grid2.deformed_element(0).vertex(0) == &grid2.deformed_element(2).vertex(3));
        REQUIRE(&grid2.deformed_element(0).vertex(1) == &grid2.deformed_element(2).vertex(1));

        grid2.add_element({1, -1});
        grid2.connect({0, 3}, {1, 0}, {0, 0});
        REQUIRE(&grid2.deformed_element(0).vertex(0) == &grid2.deformed_element(3).vertex(0));
        REQUIRE(&grid2.deformed_element(0).vertex(2) == &grid2.deformed_element(3).vertex(1));
      }
    }

    SECTION("deformed-regular")
    {
      grid3.add_element({0, -1, 0});
      grid3.add_element({0, 0, 0});
      cartdg::Regular_grid reg3 {1, 3, 0, 0.2, basis};
      reg3.add_element({1, 0, 0});
      reg3.add_element({0, 0, 1});
      grid3.connect_non_def({1, 0}, {0, 0}, {1, 0}, reg3);
      grid3.connect_non_def({1, 1}, {2, 2}, {1, 0}, reg3);
      REQUIRE(grid3.def_reg_connection(0, 0).first == &grid3.deformed_element(1));
      REQUIRE(grid3.def_reg_connection(0, 0).second == &reg3.element(0));
      REQUIRE(grid3.def_reg_connection(2, 0).first == &grid3.deformed_element(1));
      REQUIRE(grid3.def_reg_connection(2, 0).second == &reg3.element(1));
      REQUIRE_THROWS(grid3.connect_non_def({1, 0}, {1, 0}, {1, 0}, reg3));
      REQUIRE_THROWS(grid3.connect_non_def({1, 0}, {0, 0}, {0, 0}, reg3));
    }
  }

  SECTION("Jacobian calculation")
  {
    double* jac;
    grid2.add_element({0, 0});
    grid2.add_element({1, 1});
    grid2.deformed_element(0).vertex(3).pos = {0.8*0.2, 0.8*0.2, 0.};
    grid2.deformed_element(1).node_adjustments()[6 + 1] = 0.1;
    grid2.calc_jacobian();
    jac = grid2.deformed_element(0).jacobian();
    REQUIRE(jac[0*9    ] == Approx(1.));
    REQUIRE(jac[1*9    ] == Approx(0.));
    REQUIRE(jac[2*9    ] == Approx(0.));
    REQUIRE(jac[3*9    ] == Approx(1.));
    REQUIRE(jac[0*9 + 6] == Approx(1.));
    REQUIRE(jac[1*9 + 6] == Approx(-0.2));
    REQUIRE(jac[2*9 + 6] == Approx(0.));
    REQUIRE(jac[3*9 + 6] == Approx(0.8));
    REQUIRE(jac[0*9 + 8] == Approx(0.8));
    REQUIRE(jac[1*9 + 8] == Approx(-0.2));
    REQUIRE(jac[2*9 + 8] == Approx(-0.2));
    REQUIRE(jac[3*9 + 8] == Approx(0.8));
    jac = grid2.deformed_element(1).jacobian();
    REQUIRE(jac[0*9 + 5] == Approx(1.));
    REQUIRE(jac[1*9 + 5] == Approx(0.));
    REQUIRE(jac[2*9 + 5] == Approx(0.));
    REQUIRE(jac[3*9 + 5] == Approx(0.9));

    grid3.add_element({0, 0, 0});
    grid3.deformed_element(0).vertex(7).pos = {0.8*0.2, 0.8*0.2, 0.8*0.2};
    grid3.calc_jacobian();
    jac = grid3.deformed_element(0).jacobian();
    REQUIRE(jac[0] == 1.);
    REQUIRE(jac[0*27 + 26] == Approx( 0.8));
    REQUIRE(jac[1*27 + 26] == Approx(-0.2));
    REQUIRE(jac[2*27 + 26] == Approx(-0.2));
    REQUIRE(jac[7*27 + 26] == Approx(-0.2));
    REQUIRE(jac[8*27 + 26] == Approx( 0.8));
  }

  SECTION("volume integrals")
  {
    grid2.add_element({-1, 0});
    grid2.add_element({ 0, 0});
    grid2.deformed_element(1).vertex(0).pos[0] = 0.05;
    grid2.deformed_element(1).vertex(0).pos[1] = 0.07;
    grid2.calc_jacobian();
    for (int i_elem : {0, 1})
    {
      double* stage = grid2.deformed_element(i_elem).stage(0);
      for (int i_qpoint = 0; i_qpoint < grid2.n_qpoint; ++i_qpoint) stage[i_qpoint] = 1.2;
    }
    auto integral = grid2.integral();
    double area = 0.2*0.2*2. - 0.2*0.5*(0.05 + 0.07);
    REQUIRE(integral[0] == Approx(1.2*area));
  }

  SECTION("single-face integrals")
  {
    grid3.add_element({0, 0, 0});
    grid2.add_element({0, 0});
    grid2.calc_jacobian();
    grid3.calc_jacobian();
    cartdg::State_variables sv;
    SECTION("cubic polynomial, regular face")
    {
      auto pos = grid3.get_pos(0);
      double* stage = grid3.deformed_element(0).stage(0);
      for (int i_qpoint = 0; i_qpoint < grid3.n_qpoint; ++i_qpoint)
      {
        double pos1 = pos[i_qpoint + grid3.n_qpoint];
        double pos2 = pos[i_qpoint + grid3.n_qpoint*2];
        stage[i_qpoint] = std::pow(pos1, 3)*(-std::pow(pos2, 3) + pos2 + 3.) - 2.*std::pow(pos2, 2) - 1.;
      }
      REQUIRE(grid3.face_integral(sv, 0, 0, 0)[0] == Approx(-0.04081882666666667));
    }
    SECTION("constant polynomial, irregular faces")
    {
      for (int i_qpoint = 0; i_qpoint < grid3.n_qpoint; ++i_qpoint)
      {
        grid3.deformed_element(0).stage(0)[i_qpoint] = 1.;
      }
      for (int i_qpoint = 0; i_qpoint < grid2.n_qpoint; ++i_qpoint)
      {
        grid2.deformed_element(0).stage(0)[i_qpoint] = 1.;
      }

      grid3.deformed_element(0).vertex(1).pos = {0., 0., 0.8*0.2};
      grid3.calc_jacobian();
      double area = 0.2*0.2;
      REQUIRE(grid3.face_integral(sv, 0, 0, 0)[0] == Approx(0.9*area));
      REQUIRE(grid3.face_integral(sv, 0, 0, 1)[0] == Approx(area));
      REQUIRE(grid3.face_integral(sv, 0, 1, 0)[0] == Approx(0.9*area));
      REQUIRE(grid3.face_integral(sv, 0, 2, 0)[0] == Approx(area));

      grid3.deformed_element(0).vertex(3).pos = {0., 0.2, 0.8*0.2};
      grid3.calc_jacobian();
      REQUIRE(grid3.face_integral(sv, 0, 2, 1)[0] == Approx(std::sqrt(0.2*0.2 + 1.)*area));

      grid2.deformed_element(0).vertex(2).pos = {1.1*0.2, 0.};
      grid2.deformed_element(0).vertex(3).pos = {0.9*0.2, 0.9*0.2};
      grid2.deformed_element(0).vertex(0).pos = {0.3, 0.3}; // show that vertex on opposite face has no effect
      grid2.calc_jacobian();
      REQUIRE(grid2.face_integral(sv, 0, 0, 1)[0] == Approx(std::sqrt(0.2*0.2 + 0.9*0.9)*0.2));
    }
  }

  SECTION("wall surface integrals")
  {
    // some completely arbitrary elements
    grid2.add_element({ 1,0});
    grid2.add_element({-1,2});
    grid2.add_element({-1,3});
    grid2.calc_jacobian();

    grid2.add_wall(1, 0, 0);
    grid2.add_wall(2, 0, 0);
    grid2.add_wall(2, 0, 1);
    grid2.add_wall(2, 1, 1);
    for (int i_elem = 0; i_elem < grid2.n_elem; ++i_elem)
    {
      auto pos = grid2.get_pos(i_elem);
      double* state = grid2.state_r() + i_elem*grid2.n_dof;
      for (int i_qpoint = 0; i_qpoint < grid2.n_qpoint; ++i_qpoint)
      {
        state[i_qpoint] = pos[i_qpoint] + 0.1*pos[i_qpoint + grid2.n_qpoint];
      }
    }
    cartdg::State_variables sv;
    REQUIRE(grid2.surface_integral(sv)[0] == Approx(-1*.2*2*.2 - 0.2*0.2/2. + 0.1*((2*4*.2*4*.2 - 2*.2*2*.2 - 3*.2*3*.2)/2. + 4*0.2*0.2)));

    cartdg::Deformed_grid multivar (2, 3, 0, 1., basis);
    multivar.add_element({0, 0, 0});
    multivar.add_wall(0, 2, 1);
    multivar.add_wall(0, 1, 0);
    for (int i_qpoint = 0; i_qpoint < multivar.n_qpoint; ++i_qpoint)
    {
      multivar.state_r()[i_qpoint] = 2.;
      multivar.state_r()[i_qpoint + multivar.n_qpoint] = .03;
    }
    multivar.calc_jacobian();
    auto integral = multivar.surface_integral(sv);
    REQUIRE(integral[0] == Approx(4.));
    REQUIRE(integral[1] == Approx(.06));
  }

  SECTION("Adding wall boundary conditions")
  {
    grid3.add_wall(1, 2, false);
    grid3.add_wall(0, 0, true);
    REQUIRE(grid3.i_elem_wall.size() == 2);
    REQUIRE(grid3.i_dim_wall.size() == 2);
    REQUIRE(grid3.is_positive_wall.size() == 2);
    REQUIRE(grid3.i_elem_wall[0] == 1);
    REQUIRE(grid3.i_dim_wall[1] == 0);
    REQUIRE(grid3.is_positive_wall[1] == true);
  }
}
