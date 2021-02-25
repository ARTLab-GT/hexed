#include <catch.hpp>

#include <Deformed_grid.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("Vertex")
{
  cartdg::Gauss_lobatto basis (2);
  cartdg::Deformed_grid grid (1, 3, 0, 1., basis);
  grid.add_element({0, 0, 0});

  SECTION("combining")
  {
    grid.add_element({2, 2, 2});

    // setup
    cartdg::Vertex& vert1 = grid.get_vertex(1);
    REQUIRE(vert1.pos == std::array<double, 3>({0., 0., 1.}));
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 0) == 1);
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 10) == 0);
    cartdg::Vertex& vert8 = grid.get_vertex(8);
    REQUIRE(vert8.pos == std::array<double, 3>({2., 2., 2.}));
    REQUIRE(std::count(vert8.neighbor_ids.begin(), vert8.neighbor_ids.end(), 0) == 0);
    REQUIRE(std::count(vert8.neighbor_ids.begin(), vert8.neighbor_ids.end(), 10) == 1);
    REQUIRE(vert8.mass == 1);
    vert8.mass = 2; // mess just to verify that non-unit masses get added correctly
    REQUIRE(grid.vertex_ids[1] == 1);
    REQUIRE(grid.vertex_ids[8] == 8);

    // verification
    vert1.eat(vert8);
    REQUIRE(vert1.mass == 3);
    REQUIRE(vert8.mass == 0);
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 0) == 1);
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 10) == 1);
    REQUIRE(vert8.neighbor_ids.empty());
    REQUIRE(vert1.id_refs == std::vector<int>({1, 8}));
    REQUIRE(vert8.id_refs.empty());
    REQUIRE(grid.vertex_ids[1] == 1);
    REQUIRE(grid.vertex_ids[8] == 1);
    REQUIRE(vert1.pos[0] == Approx(4./3.));
    REQUIRE(vert1.pos[1] == Approx(4./3.));
    REQUIRE(vert1.pos[2] == Approx(5./3.));
  }

  SECTION("Smoothing")
  {
    cartdg::Vertex vertex = grid.get_vertex(0);
    REQUIRE(vertex.pos[0] == 0.);
    REQUIRE(vertex.pos[1] == 0.);
    REQUIRE(vertex.pos[2] == 0.);
    REQUIRE(vertex.relax[0] == 0.);
    REQUIRE(vertex.relax[1] == 0.);
    REQUIRE(vertex.relax[2] == 0.);
    vertex.calc_relax();
    REQUIRE(vertex.pos[0] == 0.);
    REQUIRE(vertex.pos[1] == 0.);
    REQUIRE(vertex.pos[2] == 0.);
    vertex.apply_relax();
    REQUIRE(vertex.pos[0] == Approx(1./6.));
    REQUIRE(vertex.pos[1] == Approx(1./6.));
    REQUIRE(vertex.pos[2] == Approx(1./6.));
    vertex.apply_relax();
    REQUIRE(vertex.pos[0] == Approx(1./6.));
    REQUIRE(vertex.pos[1] == Approx(1./6.));
    REQUIRE(vertex.pos[2] == Approx(1./6.));
    REQUIRE(vertex.relax[0] == 0.);
    REQUIRE(vertex.relax[1] == 0.);
    REQUIRE(vertex.relax[2] == 0.);
  }
}