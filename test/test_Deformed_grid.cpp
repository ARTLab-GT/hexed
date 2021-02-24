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
}