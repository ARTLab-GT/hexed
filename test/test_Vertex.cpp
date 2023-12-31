#include <catch2/catch_all.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include <hexed/Vertex.hpp>

TEST_CASE("Vertex")
{
  hexed::Vertex::Transferable_ptr ptr0 {{1.1, -3.5, 0.05}, true};
  hexed::Vertex* orig_addr = &*ptr0;
  hexed::Vertex::Transferable_ptr ptr1 {{1., 1., 1.}, false};
  REQUIRE(ptr1->is_mobile() == false);

  REQUIRE(ptr0->mass() == 1);
  REQUIRE(ptr0->pos[0] == 1.1);
  REQUIRE((*ptr0).pos[1] == -3.5);
  REQUIRE(ptr0->pos[2] == 0.05);
  REQUIRE(ptr0->is_mobile() == true);
  (*ptr0).pos[1] = 2.;
  REQUIRE(ptr0->pos[1] == 2.);
  REQUIRE(ptr0.shareable_value == 0.);
  REQUIRE(ptr1.shareable_value == 0.);
  {
    hexed::Vertex::Transferable_ptr temp(ptr0);
    REQUIRE(ptr0->mass() == 2);
  }
  REQUIRE(ptr0->mass() == 1);

  SECTION("eat")
  {
    ptr0->record.push_back(1);
    ptr0->record.push_back(-7);
    ptr1->record.push_back(4);
    ptr1->record.push_back(0);
    ptr1->record.push_back(9);
    ptr0->eat(*ptr1);
    ptr0->eat(*ptr0); // autocannibalism does nothing
    REQUIRE(ptr0->mass() == 2);
    REQUIRE(ptr0->is_mobile() == false);
    REQUIRE(ptr0->pos[0] == 1.05);
    REQUIRE(ptr0->pos[1] == 1.5);
    REQUIRE(ptr1->pos[2] == 0.525);
    ptr1->pos[0] = 2.;
    REQUIRE(&*ptr0 == orig_addr);
    REQUIRE(&*ptr1 == orig_addr);
    // records have been concatenated
    REQUIRE(ptr0->record.size() == 5);
    std::vector<double> correct {1, -7, 4, 0, 9};
    REQUIRE(std::equal(ptr0->record.begin(), ptr0->record.end(), correct.begin()));

    hexed::Vertex::Transferable_ptr ptr2 {{1., 0., 0.}, true};
    hexed::Vertex::Transferable_ptr ptr3 {{1., 0., 0.}, true};
    hexed::Vertex::Transferable_ptr ptr4 {{1., 0., 0.}, true};
    ptr2->eat(*ptr3);
    ptr3->eat(*ptr4);
    REQUIRE(ptr2->is_mobile() == true);
    ptr1->eat(*ptr3);
    REQUIRE(ptr1->mass() == 5.);
    REQUIRE(ptr0->pos[0] == 1.4);
    REQUIRE(ptr2->is_mobile() == false);
  }

  SECTION("Pointer validity")
  {
    //Note: If you choose to delete any of this test, do so with caution. Some seemingly
    //      harmless lines are capable of exposing dangling pointers in faulty implementations
    hexed::Vertex::Transferable_ptr ptr1_1 {ptr1};
    hexed::Vertex::Transferable_ptr ptr1_2 {{0., 0., 0.}};
    ptr1_2 = ptr1;

    // verifies that destructors and move semantics do not leave dangling pointers
    hexed::Vertex::Transferable_ptr* ptr1_3 = new hexed::Vertex::Transferable_ptr {ptr1};
    delete ptr1_3;
    hexed::Vertex::Transferable_ptr* ptr1_4 = new hexed::Vertex::Transferable_ptr {ptr1};
    hexed::Vertex::Transferable_ptr ptr1_5 {{0., 0., 0.}};
    ptr1_5 = std::move(*ptr1_4);
    delete ptr1_4;

    hexed::Vertex::Non_transferable_ptr ptr0_1 {*ptr0};
    hexed::Vertex::Non_transferable_ptr ptr1_6 {*ptr1};
    hexed::Vertex::Non_transferable_ptr ptr1_7 {*ptr1};
    hexed::Vertex::Non_transferable_ptr ptr1_8 {ptr1_6};
    hexed::Vertex::Non_transferable_ptr ptr1_9 {*ptr0};
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
    hexed::Vertex::Non_transferable_ptr* ptr1_10 = new hexed::Vertex::Non_transferable_ptr {*ptr1};
    delete ptr1_10;

    // explicitly nullifying works
    hexed::Vertex::Non_transferable_ptr ptr0_2 {*ptr0};
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

  SECTION("connection")
  {
    hexed::Vertex::Transferable_ptr vert0 {{0., 0., 0.}, true};
    hexed::Vertex::Transferable_ptr vert1 {{2., 0., 0.}, true};
    hexed::Vertex::Transferable_ptr vert2 {{0., 2., 0.}, true};
    hexed::Vertex::connect(*vert0, *vert1);
    hexed::Vertex::connect(*vert2, *vert0);
    // check neighbors are correct
    REQUIRE( hexed::Vertex::are_neighbors(*vert0, *vert1));
    REQUIRE( hexed::Vertex::are_neighbors(*vert1, *vert0));
    REQUIRE(!hexed::Vertex::are_neighbors(*vert1, *vert2));
    REQUIRE(!hexed::Vertex::are_neighbors(*vert2, *vert1));

    SECTION("eating and deleting")
    {
      {
        hexed::Vertex::Transferable_ptr temp {{-1., -1., -1.}, false};
        hexed::Vertex::connect(*vert0, *temp);
      }
      {
        hexed::Vertex::Transferable_ptr temp {{0., 0., 0.}, false};
        vert0->eat(*temp);
        REQUIRE(vert0->is_mobile() == false);
      }
      REQUIRE(vert0->is_mobile() == true);
      {
        hexed::Vertex::Transferable_ptr vert3 {{0., 0., 2.}, true};
        hexed::Vertex::Transferable_ptr vert4 {{-1., -1., 0.}, true};
        hexed::Vertex::connect(*vert3, *vert4);
        hexed::Vertex::connect(*vert3, *vert1);
        vert0->eat(*vert3);
        REQUIRE(vert0->mass() == 2);
      }
      REQUIRE(vert0->mass() == 1);
    }
  }

  SECTION("maximum viscosity")
  {
    ptr0->eat(*ptr1);
    hexed::Vertex::Transferable_ptr ptr2 = ptr0;
    ptr0.shareable_value = 0.1;
    ptr1.shareable_value = 0.3;
    ptr2.shareable_value = 0.2;
    REQUIRE(ptr0->shared_value() == 0.3);
    ptr0.shareable_value = 0.4;
    REQUIRE(ptr1->shared_value() == 0.4);
  }
}
