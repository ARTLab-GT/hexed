#include <catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("Deformed grid class")
{
  const int rank = std::min<int>(3, MAX_BASIS_RANK);
  cartdg::Gauss_lobatto basis (rank);
  std::vector<cartdg::Deformed_grid> grids;
  grids.emplace_back(1, 2, 0, 0.2, basis);
  grids.emplace_back(1, 3, 0, 0.2, basis);
  cartdg::Deformed_grid& grid2 = grids[0];
  cartdg::Deformed_grid& grid3 = grids[1];

  SECTION("construction")
  {
    for (cartdg::Deformed_grid grid : grids)
    {
      REQUIRE(grid.vertices.empty());
      REQUIRE(grid.vertex_ids.empty());
      REQUIRE(grid.node_adjustments.empty());
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
    grid2.add_element(std::vector<int>{1, 2});
    REQUIRE(grid2.n_elem == 2);
    REQUIRE(grid2.vertices.size() == 8);
    REQUIRE(grid2.vertex_ids.size() == 8);
    REQUIRE(grid2.node_adjustments.size() == 24);
    REQUIRE(grid2.get_vertex(0).pos[0] == Approx(-0.1));
    REQUIRE(grid2.get_vertex(0).pos[1] == Approx( 0.4));
    REQUIRE(grid2.get_vertex(1).pos[0] == Approx(-0.1));
    REQUIRE(grid2.get_vertex(1).pos[1] == Approx( 0.6));
    REQUIRE(grid2.get_vertex(7).pos[0] == Approx( 0.5));
    REQUIRE(grid2.get_vertex(7).pos[1] == Approx( 0.6));
    cartdg::Vertex& vertex0 = grid2.get_vertex(0);
    REQUIRE(vertex0.neighbor_ids.size() == 2);
    REQUIRE(std::count(vertex0.neighbor_ids.begin(), vertex0.neighbor_ids.end(), 1) == 1);
    REQUIRE(std::count(vertex0.neighbor_ids.begin(), vertex0.neighbor_ids.end(), 2) == 1);
    cartdg::Vertex& vertex5 = grid2.get_vertex(5);
    REQUIRE(vertex5.neighbor_ids.size() == 2);
    REQUIRE(std::count(vertex5.neighbor_ids.begin(), vertex5.neighbor_ids.end(), 4) == 1);
    REQUIRE(std::count(vertex5.neighbor_ids.begin(), vertex5.neighbor_ids.end(), 7) == 1);

    REQUIRE(grid3.n_elem == 0);
    grid3.add_element(std::vector<int>{2, -3, 1});
    REQUIRE(grid3.n_elem == 1);
    REQUIRE(grid3.vertices.size() == 8);
    REQUIRE(grid3.vertex_ids.size() == 8);
    REQUIRE(grid3.node_adjustments.size() == 54);
    REQUIRE(grid3.get_vertex(0).pos[0] == Approx( 0.4));
    REQUIRE(grid3.get_vertex(0).pos[1] == Approx(-0.6));
    REQUIRE(grid3.get_vertex(0).pos[2] == Approx( 0.2));
    REQUIRE(grid3.get_vertex(4).pos[0] == Approx( 0.6));
    REQUIRE(grid3.get_vertex(4).pos[1] == Approx(-0.6));
    REQUIRE(grid3.get_vertex(4).pos[2] == Approx( 0.2));
    REQUIRE(grid3.get_vertex(7).pos[0] == Approx( 0.6));
    REQUIRE(grid3.get_vertex(7).pos[1] == Approx(-0.4));
    REQUIRE(grid3.get_vertex(7).pos[2] == Approx( 0.4));
    cartdg::Vertex& vertex3 = grid3.get_vertex(3);
    REQUIRE(vertex3.neighbor_ids.size() == 3);
    REQUIRE(std::count(vertex3.neighbor_ids.begin(), vertex3.neighbor_ids.end(), 1) == 1);
    REQUIRE(std::count(vertex3.neighbor_ids.begin(), vertex3.neighbor_ids.end(), 2) == 1);
    REQUIRE(std::count(vertex3.neighbor_ids.begin(), vertex3.neighbor_ids.end(), 7) == 1);

    for (cartdg::Deformed_grid& grid : grids)
    {
      for (int i_id = 0; i_id < (int)grid.vertex_ids.size(); ++i_id)
      {
        int id = grid.vertex_ids[i_id];
        cartdg::Vertex vertex = grid.vertices[id];
        REQUIRE(vertex.parent_grid == &grid);
        REQUIRE(vertex.id == id);
        REQUIRE(vertex.id_refs.size() == 1);
        REQUIRE(vertex.id_refs[0] == i_id);
      }
    }

    SECTION("vertex interpolation")
    {
      /*
      Approximate vertex layout:
         3
                 2


      1
                   0
      */
      grid2.get_vertex(0).pos = { 0.1, 2.0, 0.0};
      grid2.get_vertex(1).pos = {-1.1, 2.1, 0.0};
      grid2.get_vertex(2).pos = { 0.0, 2.8, 0.0};
      grid2.get_vertex(3).pos = {-1.0, 2.9, 0.0};
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
      grid3.get_vertex(15).pos = {0.2*0.8, 0.2*0.8, 0.2*0.8};
      std::vector<double> pos3 = grid3.get_pos(1);
      for (int i_dim = 0; i_dim < 3; ++i_dim) REQUIRE(pos3[27*i_dim + 13] == Approx(0.2*0.475));
    }

    SECTION("warped face interpolation")
    {
      grid2.get_vertex(0).pos = {0.0, 0.0, 0.0};
      grid2.get_vertex(1).pos = {0.0, 1.0, 0.0};
      grid2.get_vertex(2).pos = {1.0, 0.0, 0.0};
      grid2.get_vertex(3).pos = {0.6, 1.0, 0.0};
      grid2.node_adjustments[2*3 + 1] =  0.2;
      grid2.node_adjustments[3*3 + 1] = -0.1;
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

      grid2.get_vertex(4).pos = {0.0, 0.0, 0.0};
      grid2.get_vertex(5).pos = {0.0, 1.0, 0.0};
      grid2.get_vertex(6).pos = {1.0, 0.0, 0.0};
      grid2.get_vertex(7).pos = {1.0, 1.0, 0.0};
      grid2.node_adjustments[4*3 + 1] =  0.1;
      std::vector<double> pos21 = grid2.get_pos(1);
      REQUIRE(pos21[ 0] == Approx(0.0));
      REQUIRE(pos21[ 6] == Approx(1.0));
      REQUIRE(pos21[ 4] == Approx(0.55));

      grid3.add_element({0, 0, 0});
      grid3.node_adjustments[6*9 + 4] = 0.01;
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
      REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[8]);
      REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[9]);
      REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[10]);
      REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[11]);

      grid3.add_element({0, 1, 0});
      grid3.connect({2, 0}, {1, 1}, {0, 1});
      REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[16]);
      REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[17]);
      REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[20]);
      REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[21]);

      grid3.add_element({0, 0, -1});
      grid3.connect({0, 3}, {2, 2}, {0, 1});
      REQUIRE(grid3.vertex_ids[0] == grid3.vertex_ids[25]);
      REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[27]);
      REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[29]);
      REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[31]);
    }
    
    SECTION("Different direction")
    {
      SECTION("0+ 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, -1, 0});
        grid3.connect({0, 1}, {0, 1}, {1, 1});
        REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[10]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[14]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[15]);
      }

      SECTION("0- 1-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({-1, 1, 0});
        grid3.connect({0, 1}, {0, 1}, {0, 0});
        REQUIRE(grid3.vertex_ids[0] == grid3.vertex_ids[8]);
        REQUIRE(grid3.vertex_ids[1] == grid3.vertex_ids[9]);
        REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[12]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[13]);
      }

      SECTION("2+ 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({0, -1, 1});
        grid3.connect({0, 1}, {2, 1}, {1, 1});
        REQUIRE(grid3.vertex_ids[1] == grid3.vertex_ids[10]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[14]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[15]);
      }

      SECTION("0+ 2+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, -1});
        grid3.connect({0, 1}, {0, 2}, {1, 1});
        REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[9]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[13]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[15]);
      }

      SECTION("2+ 0+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, -1});
        grid3.connect({1, 0}, {2, 0}, {1, 1});
        REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[9]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[13]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[15]);
      }

      SECTION("0+ 1-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({0, 1}, {0, 1}, {1, 0});
        REQUIRE(grid3.vertex_ids[4] == grid3.vertex_ids[12]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[13]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[ 8]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[ 9]);
      }

      SECTION("0- 1+")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({1, 0}, {0, 1}, {0, 1});
        REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[10]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[ 8]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[ 9]);
      }

      SECTION("1+ 0-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 1, 0});
        grid3.connect({0, 1}, {1, 0}, {1, 0});
        REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[10]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[ 8]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[ 9]);
      }

      SECTION("1+ 2-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({0, 1, 1});
        grid3.connect({0, 1}, {1, 2}, {1, 0});
        REQUIRE(grid3.vertex_ids[2] == grid3.vertex_ids[10]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[ 8]);
        REQUIRE(grid3.vertex_ids[6] == grid3.vertex_ids[14]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[12]);
      }

      SECTION("2+ 0-")
      {
        grid3.add_element({0, 0, 0});
        grid3.add_element({1, 0, 1});
        grid3.connect({0, 1}, {2, 0}, {1, 0});
        REQUIRE(grid3.vertex_ids[1] == grid3.vertex_ids[ 9]);
        REQUIRE(grid3.vertex_ids[3] == grid3.vertex_ids[11]);
        REQUIRE(grid3.vertex_ids[5] == grid3.vertex_ids[ 8]);
        REQUIRE(grid3.vertex_ids[7] == grid3.vertex_ids[10]);
        grid3.visualize("test");
      }

      SECTION("2D")
      {
        grid2.add_element({0, 0});
        grid2.add_element({1, 0});
        grid2.connect({0, 1}, {0, 0}, {1, 0});
        REQUIRE(grid2.vertex_ids[2] == grid2.vertex_ids[4]);
        REQUIRE(grid2.vertex_ids[3] == grid2.vertex_ids[5]);

        grid2.add_element({-1, -1});
        grid2.connect({0, 2}, {0, 1}, {0, 1});
        REQUIRE(grid2.vertex_ids[0] == grid2.vertex_ids[11]);
        REQUIRE(grid2.vertex_ids[1] == grid2.vertex_ids[ 9]);

        grid2.add_element({1, -1});
        grid2.connect({0, 3}, {1, 0}, {0, 0});
        REQUIRE(grid2.vertex_ids[0] == grid2.vertex_ids[12]);
        REQUIRE(grid2.vertex_ids[2] == grid2.vertex_ids[13]);
      }
    }
  }

  SECTION("Jacobian calculation")
  {
    grid2.add_element({0, 0});
    grid2.add_element({1, 1});
    grid2.get_vertex(3).pos = {0.8*0.2, 0.8*0.2, 0.};
    grid2.node_adjustments[12 + 6 + 1] = 0.1;
    grid2.calc_jacobian();
    REQUIRE(grid2.jacobian.size() == 2*4*9);
    REQUIRE(grid2.jacobian[      0*9    ] == Approx(1.));
    REQUIRE(grid2.jacobian[      1*9    ] == Approx(0.));
    REQUIRE(grid2.jacobian[      2*9    ] == Approx(0.));
    REQUIRE(grid2.jacobian[      3*9    ] == Approx(1.));
    REQUIRE(grid2.jacobian[      0*9 + 6] == Approx(1.));
    REQUIRE(grid2.jacobian[      1*9 + 6] == Approx(-0.2));
    REQUIRE(grid2.jacobian[      2*9 + 6] == Approx(0.));
    REQUIRE(grid2.jacobian[      3*9 + 6] == Approx(0.8));
    REQUIRE(grid2.jacobian[      0*9 + 8] == Approx(0.8));
    REQUIRE(grid2.jacobian[      1*9 + 8] == Approx(-0.2));
    REQUIRE(grid2.jacobian[      2*9 + 8] == Approx(-0.2));
    REQUIRE(grid2.jacobian[      3*9 + 8] == Approx(0.8));
    REQUIRE(grid2.jacobian[4*9 + 0*9 + 5] == Approx(1.));
    REQUIRE(grid2.jacobian[4*9 + 1*9 + 5] == Approx(0.));
    REQUIRE(grid2.jacobian[4*9 + 2*9 + 5] == Approx(0.));
    REQUIRE(grid2.jacobian[4*9 + 3*9 + 5] == Approx(0.9));

    grid3.add_element({0, 0, 0});
    grid3.get_vertex(7).pos = {0.8*0.2, 0.8*0.2, 0.8*0.2};
    grid3.calc_jacobian();
    REQUIRE(grid3.jacobian.size() == 9*27);
    REQUIRE(grid3.jacobian[0] == 1.);
    REQUIRE(grid3.jacobian[0*27 + 26] == Approx(0.8));
    REQUIRE(grid3.jacobian[1*27 + 26] == Approx(-0.2));
    REQUIRE(grid3.jacobian[2*27 + 26] == Approx(-0.2));
    REQUIRE(grid3.jacobian[7*27 + 26] == Approx(-0.2));
    REQUIRE(grid3.jacobian[8*27 + 26] == Approx(0.8));
  }

  SECTION("Jacobian determinant")
  {
    grid2.add_element({1, -2});
    grid2.add_element({3, 4});
    grid2.get_vertex(0).pos = {0.2, -0.39, 0.};
    grid2.get_vertex(6).pos[1] += 0.1;
    grid2.get_vertex(7).pos[1] += 0.1;
    grid2.calc_jacobian();
    REQUIRE(grid2.jacobian_det(0, 0) == Approx(0.95));
    REQUIRE(grid2.jacobian_det(0, rank - 1) == Approx(0.95));
    REQUIRE(grid2.jacobian_det(1, 0) == Approx(1.));
    REQUIRE(grid2.jacobian_det(1, rank*(rank - 1)) == Approx(1.));

    grid3.add_element({0, 0, 0});
    grid3.get_vertex(0).pos = {0.01, 0.0, 0.0};
    grid3.calc_jacobian();
    REQUIRE(grid3.jacobian_det(0, 0) == Approx(0.95));
    REQUIRE(grid3.jacobian_det(0, grid3.n_qpoint - 1) == Approx(1.));
  }

  SECTION("Integrals")
  {
    grid2.add_element({-1, 0});
    grid2.add_element({0, 0});
    grid2.get_vertex(4).pos[0] = 0.05;
    grid2.get_vertex(4).pos[1] = 0.07;
    grid2.calc_jacobian();
    for (int i_elem : {0, 1})
    {
      for (int i_qpoint = 0; i_qpoint < grid2.n_qpoint; ++i_qpoint)
      {
        grid2.state_r()[i_elem*grid2.n_dof + i_qpoint] = 1.2;
      }
    }
    auto integral = grid2.integral();
    double area = 0.2*0.2*2. - 0.2*0.5*(0.05 + 0.07);
    REQUIRE(integral[0] == Approx(1.2*area));
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