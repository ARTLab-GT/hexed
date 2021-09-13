#include <catch2/catch.hpp>

#include <Deformed_grid.hpp>
#include <Gauss_lobatto.hpp>

TEST_CASE("Vertex")
{
  cartdg::Vertex::Transferable_ptr ptr0 {{1.1, -3.5, 0.05}};
  cartdg::Vertex* orig_addr = &*ptr0;
  cartdg::Vertex::Transferable_ptr ptr1 {{1., 1., 1.}};
  ptr1->mobile = true;

  REQUIRE(ptr0->mass() == 1);
  REQUIRE(ptr0->pos[0] == 1.1);
  REQUIRE((*ptr0).pos[1] == -3.5);
  REQUIRE(ptr0->pos[2] == 0.05);
  REQUIRE(ptr0->mobile == false);
  (*ptr0).pos[1] = 2.;
  REQUIRE(ptr0->pos[1] == 2.);

  SECTION("eat")
  {
    ptr0->eat(*ptr1);
    REQUIRE(ptr0->mass() == 2);
    REQUIRE(ptr0->mobile == true);
    REQUIRE(ptr0->pos[0] == 1.05);
    REQUIRE(ptr0->pos[1] == 1.5);
    REQUIRE(ptr1->pos[2] == 0.525);
    ptr1->pos[0] = 2.;
    REQUIRE(&*ptr0 == orig_addr);
    REQUIRE(&*ptr1 == orig_addr);

    cartdg::Vertex::Transferable_ptr ptr2 {{1., 0., 0.}};
    cartdg::Vertex::Transferable_ptr ptr3 {{1., 0., 0.}};
    cartdg::Vertex::Transferable_ptr ptr4 {{1., 0., 0.}};
    ptr2->eat(*ptr3);
    ptr3->eat(*ptr4);
    ptr1->eat(*ptr3);
    REQUIRE(ptr1->mass() == 5.);
    REQUIRE(ptr0->pos[0] == 1.4);
  }

  SECTION("Pointer validity")
  {
    cartdg::Vertex::Transferable_ptr ptr1_1 {ptr1};
    cartdg::Vertex::Transferable_ptr ptr1_2 {{0., 0., 0.}};
    ptr1_2 = ptr1;

    // verifies that destructors and move semantics do not leave dangling pointers
    cartdg::Vertex::Transferable_ptr* ptr1_3 = new cartdg::Vertex::Transferable_ptr {ptr1};
    delete ptr1_3;
    cartdg::Vertex::Transferable_ptr* ptr1_4 = new cartdg::Vertex::Transferable_ptr {ptr1};
    cartdg::Vertex::Transferable_ptr ptr1_5 {{0., 0., 0.}};
    ptr1_5 = std::move(*ptr1_4);
    delete ptr1_4;

    cartdg::Vertex::Non_transferable_ptr ptr0_1 {*ptr0};
    cartdg::Vertex::Non_transferable_ptr ptr1_6 {*ptr1};
    cartdg::Vertex::Non_transferable_ptr ptr1_7 {*ptr1};
    // test that pointers are valid and point to the right Vertex
    REQUIRE(bool(ptr0_1));
    REQUIRE(&*ptr0_1 == orig_addr);
    REQUIRE(ptr0_1->pos[0] == 1.1);
    REQUIRE(bool(ptr1_6));
    REQUIRE(bool(ptr1_7));

    // destructor doesn't leave dangling pointer
    cartdg::Vertex::Non_transferable_ptr* ptr1_8 = new cartdg::Vertex::Non_transferable_ptr {*ptr1};
    delete ptr1_8;

    // explicitly nullifying works
    cartdg::Vertex::Non_transferable_ptr ptr0_2 {*ptr0};
    ptr0_2.nullify();
    REQUIRE(!ptr0_2);

    ptr0->eat(*ptr1);
    REQUIRE(&*ptr1   == orig_addr);
    REQUIRE(&*ptr1_1 == orig_addr);
    REQUIRE(&*ptr1_2 == orig_addr);
    REQUIRE(&*ptr1_5 == orig_addr);
    REQUIRE(bool(ptr0_1));
    REQUIRE(&*ptr0_1 == orig_addr);
    // `Non_transferable_ptr`s to eaten element should be nullified
    REQUIRE(!ptr1_6);
    REQUIRE(!ptr1_7);
  }
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
