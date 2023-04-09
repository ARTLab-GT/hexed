#include <catch2/catch_all.hpp>
#include <hexed/Row_index.hpp>

TEST_CASE("Row_index")
{
  hexed::Row_index ind0(3, 5, 0);
  REQUIRE(ind0.i_qpoint(0) == 0);
  REQUIRE(ind0.i_qpoint(1) == 25);
  REQUIRE(ind0.i_face_qpoint() == 0);
  ++ind0;
  REQUIRE(ind0.i_face_qpoint() == 1);
  REQUIRE(ind0.i_qpoint(0) == 1);
  REQUIRE(ind0.i_qpoint(1) == 26);
  for (int i = 0; i < 5; ++i) ++ind0;
  REQUIRE(ind0.i_face_qpoint() == 6);
  REQUIRE(ind0.i_qpoint(0) == 6);
  REQUIRE(ind0.i_qpoint(1) == 31);
  while (ind0) ++ind0;
  REQUIRE(ind0.i_face_qpoint() == 25);

  hexed::Row_index ind1(3, 5, 1);
  REQUIRE(ind1.i_qpoint(0) == 0);
  REQUIRE(ind1.i_qpoint(1) == 5);
  REQUIRE(ind1.i_face_qpoint() == 0);
  ++ind1;
  REQUIRE(ind1.i_face_qpoint() == 1);
  REQUIRE(ind1.i_qpoint(0) == 1);
  REQUIRE(ind1.i_qpoint(1) == 6);
  for (int i = 0; i < 5; ++i) ++ind1;
  REQUIRE(ind1.i_face_qpoint() == 6);
  REQUIRE(ind1.i_qpoint(0) == 26);
  REQUIRE(ind1.i_qpoint(1) == 31);
  while (ind1) ++ind1;
  REQUIRE(ind1.i_face_qpoint() == 25);

  hexed::Row_index ind2(3, 5, 2);
  REQUIRE(ind2.i_qpoint(0) == 0);
  REQUIRE(ind2.i_qpoint(1) == 1);
  REQUIRE(ind2.i_face_qpoint() == 0);
  ++ind2;
  REQUIRE(ind2.i_face_qpoint() == 1);
  REQUIRE(ind2.i_qpoint(0) == 5);
  REQUIRE(ind2.i_qpoint(1) == 6);
  for (int i = 0; i < 5; ++i) ++ind2;
  REQUIRE(ind2.i_face_qpoint() == 6);
  REQUIRE(ind2.i_qpoint(0) == 30);
  REQUIRE(ind2.i_qpoint(1) == 31);
  while (ind2) ++ind2;
  REQUIRE(ind2.i_face_qpoint() == 25);
}
