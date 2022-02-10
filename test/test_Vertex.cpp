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
  REQUIRE(ptr0.shareable_value == 0.);
  REQUIRE(ptr1.shareable_value == 0.);

  SECTION("eat")
  {
    ptr0->eat(*ptr1);
    ptr0->eat(*ptr0); // autocannibalism does nothing
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
    //Note: If you choose to delete any of this test, do so with caution. Some seemingly
    //      harmless lines are capable of exposing dangling pointers in faulty implementations
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
    cartdg::Vertex::Non_transferable_ptr ptr1_8 {ptr1_6};
    cartdg::Vertex::Non_transferable_ptr ptr1_9 {*ptr0};
    ptr1_9 = ptr1_6;
    // test that pointers are valid and point to the right Vertex
    REQUIRE(bool(ptr0_1));
    REQUIRE(&*ptr0_1 == orig_addr);
    REQUIRE(ptr0_1->pos[0] == 1.1);
    REQUIRE(bool(ptr1_6));
    REQUIRE(bool(ptr1_7));
    REQUIRE(bool(ptr1_8));
    REQUIRE(bool(ptr1_9));
    REQUIRE(&*ptr1_8 == &*ptr1_6);
    REQUIRE(&*ptr1_9 == &*ptr1_6);

    // destructor doesn't leave dangling pointer
    cartdg::Vertex::Non_transferable_ptr* ptr1_10 = new cartdg::Vertex::Non_transferable_ptr {*ptr1};
    delete ptr1_10;

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
    REQUIRE(!ptr1_8);
    REQUIRE(!ptr1_9);
  }

  SECTION("Position relaxation")
  {
    cartdg::Vertex::Transferable_ptr vert0 {{0., 0., 0.}};
    cartdg::Vertex::Transferable_ptr vert1 {{2., 0., 0.}};
    cartdg::Vertex::Transferable_ptr vert2 {{0., 2., 0.}};
    cartdg::Vertex::connect(*vert0, *vert1);
    cartdg::Vertex::connect(*vert2, *vert0);

    SECTION("calling once")
    {
      vert0->calc_relax();
      REQUIRE(vert0->pos == std::array<double, 3>{0., 0., 0.});
      vert0->apply_relax();
      REQUIRE(vert0->pos == std::array<double, 3>{0., 0., 0.});
      vert0->mobile = true;
      vert0->calc_relax();
      vert0->apply_relax();
      REQUIRE(vert1->pos == std::array<double, 3>{2., 0., 0.});
      REQUIRE(vert0->pos == std::array<double, 3>{.5, .5, 0.});
    }
    vert0->mobile = true;

    SECTION("calling multiple times")
    {
      vert0->calc_relax();
      vert0->calc_relax();
      vert0->apply_relax();
      vert0->apply_relax();
      REQUIRE(vert1->pos == std::array<double, 3>{2., 0., 0.});
      REQUIRE(vert0->pos == std::array<double, 3>{.5, .5, 0.});
    }

    SECTION("eating and deleting")
    {
      cartdg::Vertex::Transferable_ptr* temp = new cartdg::Vertex::Transferable_ptr {{-1., -1., -1.}};
      cartdg::Vertex::connect(*vert0, **temp);
      delete temp;
      cartdg::Vertex::Transferable_ptr vert3 {{0., 0., 2.}};
      cartdg::Vertex::Transferable_ptr vert4 {{-1., -1., 0.}};
      cartdg::Vertex::connect(*vert3, *vert4);
      cartdg::Vertex::connect(*vert3, *vert1);
      vert0->eat(*vert3);
      vert0->calc_relax();
      vert0->apply_relax();
      REQUIRE(vert0->pos == std::array<double, 3>{1./6., 1./6., .5});
    }

    SECTION("self-connection does nothing")
    {
      cartdg::Vertex::connect(*vert0, *vert0);
      vert0->calc_relax();
      vert0->apply_relax();
      REQUIRE(vert0->pos == std::array<double, 3>{.5, .5, 0.});
    }
  }

  SECTION("maximum viscosity")
  {
    ptr0->eat(*ptr1);
    cartdg::Vertex::Transferable_ptr ptr2 = ptr0;
    ptr0.shareable_value = 0.1;
    ptr1.shareable_value = 0.3;
    ptr2.shareable_value = 0.2;
    REQUIRE(ptr0->shared_max_value() == 0.3);
    ptr0.shareable_value = 0.4;
    REQUIRE(ptr1->shared_max_value() == 0.4);
    ptr0.shareable_value = ptr1.shareable_value = ptr2.shareable_value = -0.1;
    REQUIRE(ptr0->shared_max_value() == 0.);
  }
}
