#include <catch2/catch.hpp>

#include <Deformed_grid.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("Vertex")
{
  cartdg::Vertex::Transferable_ptr ptr0 {{1.1, -3.5, 0.05}};
  REQUIRE(ptr0->mass() == 1);
  REQUIRE(ptr0->pos[0] == 1.1);
  REQUIRE((*ptr0).pos[1] == -3.5);
  REQUIRE(ptr0->pos[2] == 0.05);
  REQUIRE(ptr0->mobile == false);
  (*ptr0).pos[1] = 2.;
  REQUIRE(ptr0->pos[1] == 2.);

  cartdg::Vertex::Transferable_ptr ptr1 {{1., 1., 1.}};
  ptr0->eat(*ptr1);
  REQUIRE(ptr0->mass() == 2);
  REQUIRE(ptr1->mass() == 2);
  REQUIRE(ptr0->pos[0] == 1.05);
  REQUIRE(ptr0->pos[1] == 1.5);
  REQUIRE(ptr1->pos[2] == 0.525);
  ptr0->pos[0] = 4.;
  REQUIRE(ptr1->pos[0] == 4.);
  ptr1->pos[0] = 2.;
  REQUIRE(ptr0->pos[0] == 2.);

  cartdg::Vertex::Transferable_ptr ptr2 {{1., 0., 0.}};
  cartdg::Vertex::Transferable_ptr ptr3 {{1., 0., 0.}};
  cartdg::Vertex::Transferable_ptr ptr4 {{1., 0., 0.}};
  ptr2->eat(*ptr3);
  ptr3->eat(*ptr4);
  ptr1->eat(*ptr3);
  REQUIRE(ptr1->mass() == 5.);
  REQUIRE(ptr0->pos[0] == 1.4);
}

TEST_CASE("Vertex (old interface)")
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
    vert8.mobile = true;
    REQUIRE(vert8.pos == std::array<double, 3>({2., 2., 2.}));
    REQUIRE(std::count(vert8.neighbor_ids.begin(), vert8.neighbor_ids.end(), 0) == 0);
    REQUIRE(std::count(vert8.neighbor_ids.begin(), vert8.neighbor_ids.end(), 10) == 1);
    REQUIRE(vert8.mass() == 1);
    REQUIRE(grid.vertex_ids[1] == 1);
    REQUIRE(grid.vertex_ids[8] == 8);
    REQUIRE(vert1.mobile == false);

    // verification
    cartdg::Vertex another8(1000);
    another8.pos = vert8.pos;
    another8.parent_grid = vert8.parent_grid;
    vert8.eat(another8); // want to make sure following works when vert8.mass() != 1
    vert1.eat(vert8);
    REQUIRE(vert1.mass() == 3);
    REQUIRE(vert8.mass() == 0);
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
    REQUIRE(vert1.mobile == true);

    vert1.eat(vert1);
    REQUIRE(vert1.mass() == 3);
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 0) == 1);
    REQUIRE(std::count(vert1.neighbor_ids.begin(), vert1.neighbor_ids.end(), 10) == 1);
    REQUIRE(vert1.id_refs == std::vector<int>({1, 8}));
    REQUIRE(grid.vertex_ids[1] == 1);
    REQUIRE(vert1.pos[0] == Approx(4./3.));
    REQUIRE(vert1.pos[1] == Approx(4./3.));
    REQUIRE(vert1.pos[2] == Approx(5./3.));
    REQUIRE(vert1.mobile == true);

    cartdg::Deformed_grid other_grid (1, 3, 0, 1., basis);
    other_grid.add_element({0, 0, 0});
    REQUIRE_THROWS(grid.get_vertex(0).eat(other_grid.get_vertex(0)));
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
    REQUIRE(vertex.pos[0] == 0.);
    REQUIRE(vertex.pos[1] == 0.);
    REQUIRE(vertex.pos[2] == 0.);
    REQUIRE(vertex.relax[0] == 0.);
    REQUIRE(vertex.relax[1] == 0.);
    REQUIRE(vertex.relax[2] == 0.);

    vertex.mobile = true;
    vertex.calc_relax();
    REQUIRE(vertex.pos[0] == 0.);
    REQUIRE(vertex.pos[1] == 0.);
    REQUIRE(vertex.pos[2] == 0.);
    vertex.apply_relax();
    REQUIRE(vertex.pos[0] == Approx(1./6.));
    REQUIRE(vertex.pos[1] == Approx(1./6.));
    REQUIRE(vertex.pos[2] == Approx(1./6.));
    REQUIRE(vertex.relax[0] == 0.);
    REQUIRE(vertex.relax[1] == 0.);
    REQUIRE(vertex.relax[2] == 0.);
    vertex.apply_relax();
    REQUIRE(vertex.pos[0] == Approx(1./6.));
    REQUIRE(vertex.pos[1] == Approx(1./6.));
    REQUIRE(vertex.pos[2] == Approx(1./6.));
  }
}
