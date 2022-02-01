#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Regular_grid.hpp>
#include <Gauss_lobatto.hpp>
#include <Gauss_legendre.hpp>
#include <math.hpp>
#include <Tecplot_file.hpp>

TEST_CASE("Deformed grid class")
{
  const int row_size = std::min<int>(3, cartdg::config::max_row_size);
  cartdg::Gauss_lobatto basis (row_size);
  std::vector<cartdg::Deformed_grid> grids;
  grids.emplace_back(1, 2, 0, 0.2, basis);
  grids.emplace_back(1, 3, 0, 0.2, basis);
  cartdg::Deformed_grid& grid2 = grids[0];
  cartdg::Deformed_grid& grid3 = grids[1];
  REQUIRE_THROWS(cartdg::Deformed_grid (1, 2, 1, 0.2, basis));

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
      cartdg::Gauss_legendre leg_basis {row_size};
      cartdg::Deformed_grid grid {1, 2, 0, 0.2, leg_basis};
      grid.origin[0] = 0.1;
      grid.add_element({-1, 2});
      /*
      Approximate vertex layout:
         3
                 2


      1
                   0
      */
      grid.deformed_element(0).vertex(0).pos = { 0.1, 2.0, 0.0};
      grid.deformed_element(0).vertex(1).pos = {-1.1, 2.1, 0.0};
      grid.deformed_element(0).vertex(2).pos = { 0.0, 2.8, 0.0};
      grid.deformed_element(0).vertex(3).pos = {-1.0, 2.9, 0.0};
      std::vector<double> pos2 = grid.get_pos(0);
      Eigen::MatrixXd bound_mat = grid.basis.boundary();
      std::vector<double> bound_pos (8);
      for (int i_dim : {0, 1})
      {
        const int n_qpoint {row_size*row_size};
        Eigen::Map<Eigen::VectorXd> bound_pos_mat {bound_pos.data() + 4*i_dim, 4};
        Eigen::Map<Eigen::VectorXd> qpoint_pos_mat {pos2.data() + n_qpoint*i_dim, n_qpoint};
        bound_pos_mat = cartdg::custom_math::hypercube_matvec(bound_mat, qpoint_pos_mat);
      }
      REQUIRE(bound_pos[0] == Approx( 0.1));
      REQUIRE(bound_pos[1] == Approx(-1.1));
      REQUIRE(bound_pos[4] == Approx( 2. ));
      REQUIRE(bound_pos[6] == Approx( 2.8));
      REQUIRE(bound_pos[7] == Approx( 2.9));

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

      CHECK(pos2[ 3] == Approx(0.5 - 0.2*0.2));
      CHECK(pos2[ 4] == Approx(0.4 - 0.2*(0.2 - 0.1)/2));
      CHECK(pos2[ 5] == Approx(0.3 + 0.2*0.1));
      CHECK(pos2[12] == Approx(0.0 + 0.2));
      CHECK(pos2[13] == Approx(0.5 + (0.2 - 0.1)/2));
      CHECK(pos2[14] == Approx(1.0 - 0.1));

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

      cartdg::Gauss_legendre leg_basis {row_size};
      cartdg::Deformed_grid leg_grid {1, 2, 0, 0.2, leg_basis};
      leg_grid.add_element({0, 0});
      leg_grid.deformed_element(0).node_adjustments()[1] = 0.1;
      leg_grid.deformed_element(0).node_adjustments()[3] = -0.2;
      std::vector<double> pos {leg_grid.get_pos(0)};
      REQUIRE(pos[3] == Approx(0.08));
      REQUIRE(pos[4] == Approx(0.11));
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

      REQUIRE(grid3.connection(1).face_index(0).element == &grid3.deformed_element(2));
      REQUIRE(grid3.connection(1).face_index(1).element == &grid3.deformed_element(0));
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

        REQUIRE(grid3.connection(0).face_index(0).i_dim == 0);
        REQUIRE(grid3.connection(0).face_index(1).i_dim == 1);
        REQUIRE(grid3.connection(0).face_index(0).is_positive == true);
        REQUIRE(grid3.connection(0).face_index(1).is_positive == false);
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

    SECTION("interface synch")
    {
      cartdg::Gauss_lobatto lin_basis {2};
      SECTION("3D")
      {
        cartdg::Deformed_grid lin_grid {5, 3, 0, 1., lin_basis};
        lin_grid.add_element({0, 0, 0});
        lin_grid.add_element({1, 0, 0});
        lin_grid.connect({0, 1}, {0, 0}, {1, 0});
        {
          // make the qpoints align imperfectly to check that the Jacobian is properly averaged
          auto& elem = lin_grid.deformed_element(1);
          elem.node_adjustments()[(2*2 + 0)*elem.storage_params().n_qpoint()/2 + 0] = 0.3;
        }
        // connect different dimensions to verify qpoint ordering and axis permutation
        lin_grid.add_element({0, 0, 0});
        lin_grid.add_element({1, 1, 0});
        lin_grid.connect({2, 3}, {0, 1}, {1, 0});
        lin_grid.add_element({1, -1, 0});
        lin_grid.add_element({0, 0, 0});
        lin_grid.connect({4, 5}, {1, 0}, {1, 1});
        lin_grid.add_element({0, 0, 0});
        lin_grid.add_element({-1, 0, -1});
        lin_grid.connect({6, 7}, {0, 2}, {0, 1});
        lin_grid.calc_jacobian();
        REQUIRE(lin_grid.connection(0).jacobian(0, 0, 0) == Approx(1.));
        REQUIRE(lin_grid.connection(0).jacobian(2, 1, 0) == Approx(-0.3/2.));
        REQUIRE(lin_grid.connection(0).jacobian(2, 2, 0) == Approx(1. - 0.3/2.));
        REQUIRE(lin_grid.connection(1).jacobian(0, 0, 0) == Approx(1.));
        REQUIRE(lin_grid.connection(1).jacobian(0, 0, 2) == Approx(0.5));
        REQUIRE(lin_grid.connection(1).jacobian(0, 1, 2) == Approx(-0.5));
        REQUIRE(lin_grid.connection(1).jacobian(1, 1, 2) == Approx(0.5));
        REQUIRE(lin_grid.connection(2).jacobian(0, 0, 0) == Approx(0.5));
        REQUIRE(lin_grid.connection(2).jacobian(0, 1, 0) == Approx(-0.5));
        REQUIRE(lin_grid.connection(2).jacobian(2, 2, 0) == Approx(1.));
        REQUIRE(lin_grid.connection(3).jacobian(0, 0, 0) == Approx(-0.5));
        REQUIRE(lin_grid.connection(3).jacobian(0, 0, 1) == Approx(-0.75));
      }
      SECTION("2D")
      {
        cartdg::Deformed_grid lin_grid {4, 2, 0, 1., lin_basis};
        lin_grid.add_element({0, 0, 0});
        lin_grid.add_element({1, 1, 0});
        lin_grid.connect({0, 1}, {0, 1}, {1, 0});
        {
          // make the qpoints align imperfectly to check that the Jacobian is properly averaged
          auto& elem = lin_grid.deformed_element(0);
          elem.node_adjustments()[(1*2 + 1)*elem.storage_params().n_qpoint()/2 + 1] = 0.1;
        }
        // connect different dimensions to verify qpoint ordering and axis permutation
        lin_grid.calc_jacobian();
        REQUIRE(lin_grid.connection(0).jacobian(0, 0, 0) == Approx(1.));
        REQUIRE(lin_grid.connection(0).jacobian(0, 1, 0) == Approx(-0.5 - 0.05/std::sqrt(2.)));
        REQUIRE(lin_grid.connection(0).jacobian(0, 0, 1) == Approx(0.5 - 0.05*std::sqrt(2.)));
      }
    }
  }

  SECTION("vertex relaxation")
  {
    grid3.add_element({0, 0, 0});
    grid3.add_element({1, 0, 0});
    grid3.connect({0, 1}, {0, 0}, {1, 0});
    for (int i_elem : {0, 1})
    {
      for (int i_vertex = 0; i_vertex < 8; ++i_vertex)
      {
        grid3.deformed_element(i_elem).vertex(i_vertex).mobile = true;
      }
    }
    SECTION("without purging")
    {
      grid3.calc_vertex_relaxation();
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[0] == 0.);
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[1] == 0.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      grid3.apply_vertex_relaxation();
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[0] == 0.2*0.5/3.);
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[1] == 0.2*0.5/3.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[1] == 0.2*0.5/4.);
    }
    SECTION("with purging")
    {
      grid3.purge_vertices();
      grid3.calc_vertex_relaxation();
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[0] == 0.);
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[1] == 0.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      grid3.apply_vertex_relaxation();
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[0] == 0.2*0.5/3.);
      REQUIRE(grid3.deformed_element(0).vertex(0).pos[1] == 0.2*0.5/3.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[0] == 0.2*1.);
      REQUIRE(grid3.deformed_element(0).vertex(4).pos[1] == 0.2*0.5/4.);
    }
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
    cartdg::Surface_from_domain sfd {sv};
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
      REQUIRE(grid3.face_integral(sfd, 0, 0, 0)[0] == Approx(-0.04081882666666667));
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
      REQUIRE(grid3.face_integral(sfd, 0, 0, 0)[0] == Approx(0.9*area));
      REQUIRE(grid3.face_integral(sfd, 0, 0, 1)[0] == Approx(area));
      REQUIRE(grid3.face_integral(sfd, 0, 1, 0)[0] == Approx(0.9*area));
      REQUIRE(grid3.face_integral(sfd, 0, 2, 0)[0] == Approx(area));

      grid3.deformed_element(0).vertex(3).pos = {0., 0.2, 0.8*0.2};
      grid3.calc_jacobian();
      REQUIRE(grid3.face_integral(sfd, 0, 2, 1)[0] == Approx(std::sqrt(0.2*0.2 + 1.)*area));

      grid2.deformed_element(0).vertex(2).pos = {1.1*0.2, 0.};
      grid2.deformed_element(0).vertex(3).pos = {0.9*0.2, 0.9*0.2};
      grid2.deformed_element(0).vertex(0).pos = {0.3, 0.3}; // show that vertex on opposite face has no effect
      grid2.calc_jacobian();
      REQUIRE(grid2.face_integral(sfd, 0, 0, 1)[0] == Approx(std::sqrt(0.2*0.2 + 0.9*0.9)*0.2));
    }
  }

  SECTION("wall surface integrals")
  {
    cartdg::Gauss_legendre legendre {row_size}; // make sure that it works with legendre basis
    cartdg::Deformed_grid grid {1, 2, 0, 0.2, legendre};

    // some completely arbitrary elements
    grid.add_element({ 1,0});
    grid.add_element({-1,2});
    grid.add_element({-1,3});
    grid.calc_jacobian();

    grid.add_wall(1, 0, 0);
    grid.add_wall(2, 0, 0);
    grid.add_wall(2, 0, 1);
    grid.add_wall(2, 1, 1);
    REQUIRE(grid.def_elem_wall(1).face_index().element == &grid.deformed_element(2));
    REQUIRE(grid.def_elem_wall(1).i_elem() == 2);
    REQUIRE(grid.def_elem_wall(1).face_index().i_dim == 0);
    REQUIRE(grid.def_elem_wall(3).face_index().i_dim == 1);
    REQUIRE(grid.def_elem_wall(2).face_index().is_positive == true);

    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      auto pos = grid.get_pos(i_elem);
      double* stage = grid.deformed_element(i_elem).stage(0);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        stage[i_qpoint] = pos[i_qpoint] + 0.1*pos[i_qpoint + grid.n_qpoint];
      }
    }
    cartdg::State_variables sv;
    cartdg::Surface_from_domain sfd {sv};
    REQUIRE(grid.surface_integral(sfd)[0] == Approx(-1*.2*2*.2 - 0.2*0.2/2. + 0.1*((2*4*.2*4*.2 - 2*.2*2*.2 - 3*.2*3*.2)/2. + 4*0.2*0.2)));

    cartdg::Deformed_grid multivar (2, 3, 0, 1., basis);
    multivar.add_element({0, 0, 0});
    multivar.add_wall(0, 2, 1);
    multivar.add_wall(0, 1, 0);
    for (int i_qpoint = 0; i_qpoint < multivar.n_qpoint; ++i_qpoint)
    {
      multivar.element(0).stage(0)[i_qpoint] = 2.;
      multivar.element(0).stage(0)[i_qpoint + multivar.n_qpoint] = .03;
    }
    multivar.calc_jacobian();
    auto integral = multivar.surface_integral(sfd);
    REQUIRE(integral[0] == Approx(4.));
    REQUIRE(integral[1] == Approx(.06));

    cartdg::Deformed_grid warped {4, 2, 0, 1., legendre};
    warped.add_element({0, 0});
    warped.deformed_element(0).vertex(0).pos[0] = -2.;
    warped.calc_jacobian();
    warped.add_wall(0, 0, 0);
    warped.add_wall(0, 1, 1);
    double* stage {warped.element(0).stage(0)};
    for (int i_qpoint = 0; i_qpoint < warped.n_qpoint; ++i_qpoint) {
      stage[i_qpoint + 0*warped.n_qpoint] = 0.;
      stage[i_qpoint + 1*warped.n_qpoint] = 0.;
      stage[i_qpoint + 2*warped.n_qpoint] = 1.3;
      stage[i_qpoint + 3*warped.n_qpoint] = 1e5/0.4;
    }
    cartdg::Force_per_area fpa;
    auto force_integral {warped.surface_integral(fpa)};
    REQUIRE(force_integral[0] == Approx(-1e5));
    REQUIRE(force_integral[1] == Approx( 3e5));
  }

  SECTION("degenerate handling")
  {
    int row_size = cartdg::config::max_row_size;
    cartdg::Gauss_legendre gleg (row_size);
    cartdg::Deformed_grid grid (1, 2, 0, 1., gleg);
    for (int i_elem = 0; i_elem < 2; ++i_elem) {
      grid.add_element({i_elem, 0, 0});
      grid.deformed_element(i_elem).vertex(2).pos[0] -= 1.;
      grid.deformed_element(i_elem).vertex(3).pos[1] -= 0.1;
      grid.deformed_element(i_elem).vertex(1).pos[0] += 0.3;
    }
    grid.calc_jacobian();
    grid.deformed_element(1).degenerate = true;
    srand(10);
    SECTION("random") // convenient visualization of projection
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
        grid.element(0).stage(0)[i_qpoint] = grid.element(1).stage(0)[i_qpoint] = (rand()%1000)/10000.;
      }
      double integral = grid.integral()[0];
      grid.project_degenerate(0);
      cartdg::Tecplot_file file {"projection_random", 2, 1, 0.};
      grid.visualize_interior(file);
      REQUIRE(grid.integral()[0] - integral == Approx(0.).scale(1.)); // test conservation
    }
    SECTION("perturbation") // make sure information near degenerate point is annihilated
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
        // use stage 2 just to make sure you can
        grid.element(0).stage(2)[i_qpoint] = grid.element(1).stage(2)[i_qpoint] = 0.;
      }
      grid.element(0).stage(2)[row_size] = grid.element(1).stage(2)[row_size] = 1.;
      grid.element(0).stage(2)[(row_size - 2)*row_size] = grid.element(1).stage(2)[(row_size - 2)*row_size] = -1.;
      grid.project_degenerate(2);
      // compute a sort of uniform norm on the quadrature points to assess the amount of
      // oscillation suppression
      double elem_norm [2] {};
      for (int i_elem = 0; i_elem < 2; ++i_elem) {
        for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
          elem_norm[i_elem] = std::max(elem_norm[i_elem], std::abs(grid.element(i_elem).stage(2)[i_qpoint]));
          grid.element(i_elem).stage(0)[i_qpoint] = grid.element(i_elem).stage(2)[i_qpoint]; // for plotting
        }
      }
      cartdg::Tecplot_file file {"projection_perturbation", 2, 1, 0.};
      grid.visualize_interior(file);
      REQUIRE(elem_norm[0]/elem_norm[1] > 10.);
    }
    SECTION("gaussian")
    {
      std::vector<double> diffs;
      for (int resolution : {1, 2}) {
        for (int i_elem = 0; i_elem < 2; ++i_elem) {
          std::vector<double> pos = grid.get_pos(i_elem);
          for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
            double x0 = (pos[i_qpoint] - i_elem)/resolution;
            double x1 = pos[grid.n_qpoint + i_qpoint]/resolution;
            grid.element(i_elem).stage(0)[i_qpoint] = std::exp(-(x0*x0 + x1*x1));
          }
        }
        grid.project_degenerate(0);
        if (resolution == 1) {
          cartdg::Tecplot_file file {"projection_gaussian", 2, 1, 0.};
          grid.visualize_interior(file);
        }
        // compute the L2 norm of the difference with/without projection. This is a little
        // awkward because *someone* wrote an integral function that can't do individual elements
        for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
          double diff = grid.element(0).stage(0)[i_qpoint] - grid.element(1).stage(0)[i_qpoint];
          grid.element(0).stage(0)[i_qpoint] = diff*diff;
          grid.element(1).stage(0)[i_qpoint] = 0.;
        }
        diffs.push_back(std::sqrt(grid.integral()[0]));
      }
      // demand a greater order of accuracy than required
      REQUIRE(diffs[0]/diffs[1] > std::pow(2, row_size));
    }
  }
}
