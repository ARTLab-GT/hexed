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
  SECTION("add_element")
  {
    REQUIRE(grid2.n_elem == 0);
    grid2.origin[0] = 0.1;
    grid2.add_element(std::vector<int>{-1, 2});
    grid2.add_element(std::vector<int>{1, 2});
    REQUIRE(grid2.n_elem == 2);
    REQUIRE(grid2.vertices.size() == 8);
    REQUIRE(grid2.vertex_ids.size() == 8);
    REQUIRE(grid2.node_adjustments.size() == 24);
    REQUIRE(grid2.vertices[grid2.vertex_ids[0]].pos[0] == Approx(-0.1));
    REQUIRE(grid2.vertices[grid2.vertex_ids[0]].pos[1] == Approx( 0.4));
    REQUIRE(grid2.vertices[grid2.vertex_ids[1]].pos[0] == Approx(-0.1));
    REQUIRE(grid2.vertices[grid2.vertex_ids[1]].pos[1] == Approx( 0.6));
    REQUIRE(grid2.vertices[grid2.vertex_ids[7]].pos[0] == Approx( 0.5));
    REQUIRE(grid2.vertices[grid2.vertex_ids[7]].pos[1] == Approx( 0.6));
    cartdg::Vertex& vertex0 = grid2.vertices[grid2.vertex_ids[0]];
    REQUIRE(vertex0.neighbor_ids.size() == 2);
    REQUIRE(std::count(vertex0.neighbor_ids.begin(), vertex0.neighbor_ids.end(), 1) == 1);
    REQUIRE(std::count(vertex0.neighbor_ids.begin(), vertex0.neighbor_ids.end(), 2) == 1);
    cartdg::Vertex& vertex5 = grid2.vertices[grid2.vertex_ids[5]];
    REQUIRE(vertex5.neighbor_ids.size() == 2);
    REQUIRE(std::count(vertex5.neighbor_ids.begin(), vertex5.neighbor_ids.end(), 4) == 1);
    REQUIRE(std::count(vertex5.neighbor_ids.begin(), vertex5.neighbor_ids.end(), 7) == 1);

    REQUIRE(grid3.n_elem == 0);
    grid3.add_element(std::vector<int>{2, -3, 1});
    REQUIRE(grid3.n_elem == 1);
    REQUIRE(grid3.vertices.size() == 8);
    REQUIRE(grid3.vertex_ids.size() == 8);
    REQUIRE(grid3.node_adjustments.size() == 54);
    REQUIRE(grid3.vertices[grid3.vertex_ids[0]].pos[0] == Approx( 0.4));
    REQUIRE(grid3.vertices[grid3.vertex_ids[0]].pos[1] == Approx(-0.6));
    REQUIRE(grid3.vertices[grid3.vertex_ids[0]].pos[2] == Approx( 0.2));
    REQUIRE(grid3.vertices[grid3.vertex_ids[4]].pos[0] == Approx( 0.6));
    REQUIRE(grid3.vertices[grid3.vertex_ids[4]].pos[1] == Approx(-0.6));
    REQUIRE(grid3.vertices[grid3.vertex_ids[4]].pos[2] == Approx( 0.2));
    REQUIRE(grid3.vertices[grid3.vertex_ids[7]].pos[0] == Approx( 0.6));
    REQUIRE(grid3.vertices[grid3.vertex_ids[7]].pos[1] == Approx(-0.4));
    REQUIRE(grid3.vertices[grid3.vertex_ids[7]].pos[2] == Approx( 0.4));
    cartdg::Vertex& vertex3 = grid3.vertices[grid2.vertex_ids[3]];
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
  }
}