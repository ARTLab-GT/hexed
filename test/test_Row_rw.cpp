#include <catch2/catch.hpp>
#include <hexed/Row_rw.hpp>

TEST_CASE("Row_rw")
{
  double interior [3*7*7];
  for (int i = 0; i < 3*7*7; ++i) interior[i] = i;
  hexed::Row_index ind(2, 7, 0);
  ++ind; ++ind;
  auto row = hexed::Row_rw<3, 7>::read_row(interior, ind);
  REQUIRE(row.rows() == 7);
  REQUIRE(row.cols() == 3);
  REQUIRE(row(0) == Approx(2));
  REQUIRE(row(6) == Approx(44));
  REQUIRE(row(7) == Approx(51));
  ++ind;
  hexed::Row_rw<3, 7>::write_row(row, interior, ind, .1);
  REQUIRE(interior[3] == Approx(2.3));
  REQUIRE(interior[52] == Approx(51 + 5.2));
  REQUIRE(interior[17] == Approx(16 + 1.7));
  REQUIRE(interior[2] == Approx(2));

  double face[4*7*2*2];
  for (int i_face = 0; i_face < 4; ++i_face) {
    for (int i = 0; i < 4*7; ++i) {
      face[i_face*4*7 + i] = i + .1*i_face;
    }
  }
  auto bound = hexed::Row_rw<3, 7>::read_bound(face, ind);
  REQUIRE(bound.rows() == 2);
  REQUIRE(bound.cols() == 3);
  REQUIRE(bound(0, 0) == Approx(3.0));
  REQUIRE(bound(1, 0) == Approx(3.1));
  REQUIRE(bound(0, 1) == Approx(10.0));
  bound = hexed::Row_rw<3, 7>::read_bound(face, hexed::Row_index(2, 7, 1));
  REQUIRE(bound(0, 0) == Approx(0.2));
  REQUIRE(bound(1, 0) == Approx(0.3));
  hexed::Row_rw<3, 7>::write_bound(Eigen::Matrix<double, 2, 3>::Ones(), face, ind);
  REQUIRE(face[ 3] == Approx(1));
  REQUIRE(face[31] == Approx(1));
}
