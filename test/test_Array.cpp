#include <catch2/catch_all.hpp>
#define HEXED_ARRAY_BOUNDS_CHECK true
#include <hexed/Array.hpp>

hexed::Array<double> make_array()
{
  hexed::Array<double> arr({2});
  arr[0] = .0;
  arr[1] = .1;
  return arr;
}

TEST_CASE("Array")
{
  hexed::Array<double> arr0({2, 3, 4});
  REQUIRE(arr0.order() == 3);
  REQUIRE_THAT(arr0.shape(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3, 4}));
  REQUIRE(arr0.size() == 24);

  hexed::Array<double> zero_size({});
  REQUIRE(zero_size.shape().empty());
  REQUIRE(zero_size.size() == 0);
  REQUIRE(arr0.data()[0] == 0.);
  REQUIRE(arr0[23] == 0.);

  auto view = arr0();
  view[1] = 0.6;
  REQUIRE(arr0[1] == Catch::Approx(0.6));
  arr0[20] = -1.;
  REQUIRE(arr0(0).order() == 2);
  REQUIRE(arr0(1)(0).size() == 4);
  REQUIRE_THAT(arr0(1).shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 4}));
  REQUIRE(arr0(1)(2)[0] == Catch::Approx(-1.));
  REQUIRE_THROWS(arr0(0)(0)(0)(0));
  REQUIRE_THROWS(arr0(2));
  REQUIRE_THROWS(arr0(0)(0)(0)[0]);
  REQUIRE_THROWS(arr0[24]);
  hexed::Array<double> arr1({2, 3, 4});
  for (int i = 0; i < 24; ++i) arr1[i] = i;
  arr0 = arr1;
  for (int i = 0; i < 24; ++i) REQUIRE(arr0[i] == Catch::Approx(i));
  arr0[1] = 42;
  REQUIRE(arr1[1] == Catch::Approx(1));

  hexed::Array<double> arr2(arr0);
  REQUIRE_THAT(arr0.shape(), Catch::Matchers::RangeEquals(std::vector<int>{2, 3, 4}));
  REQUIRE(arr2[1] == Catch::Approx(42));
  arr2[1] = 406;
  REQUIRE(arr0[1] == Catch::Approx(406));
  hexed::Array<double> arr3(std::move(make_array()));
  REQUIRE(arr3[1] == Catch::Approx(.1));
  hexed::Array<double> arr4(arr0.copy());
  arr4[1] = 287.0528;
  REQUIRE(arr0[1] == Catch::Approx(406));

  REQUIRE(arr0.same_shape(arr1));
  REQUIRE(arr1.same_shape(arr0));
  hexed::Array<double> arr5({3, 4, 2});
  REQUIRE(!arr0.same_shape(arr5));
  REQUIRE(!arr5.same_shape(arr0));

  hexed::Array<double> arr6({5, 2});
  hexed::Array<double> arr7(arr6(1, 4));
  REQUIRE_THAT(arr7.shape(), Catch::Matchers::RangeEquals(std::vector<int>{3, 2}));
  REQUIRE(&arr7(0)[0] == &arr6(1)[0]);
  REQUIRE(&arr7(2)[1] == &arr6(3)[1]);
  REQUIRE_THAT(arr6(1, 10).shape(), Catch::Matchers::RangeEquals(std::vector<int>{4, 2}));
  REQUIRE_THAT(arr6(1,  0).shape(), Catch::Matchers::RangeEquals(std::vector<int>{0, 2}));
  REQUIRE_THAT(arr6(6,  7).shape(), Catch::Matchers::RangeEquals(std::vector<int>{0, 2}));
  REQUIRE_THROWS(arr6(0)(0)(0, 1));

  auto arr8{hexed::Array<double>::make(.1, -.2, 1.5)};
  REQUIRE_THAT(arr8, Catch::Matchers::RangeEquals(std::vector<double>{.1, -.2, 1.5}, hexed::math::Approx_equal()));

  SECTION("arithmetic")
  {
    // only test one binary operator, since all of them are defined with the same macro
    hexed::Array<double> a0({2, 2});
    a0[0] = .0;
    a0[1] = .1;
    a0[2] = .2;
    a0[3] = .3;
    hexed::Array<double> a1({4});
    a1[0] = 0.;
    a1[1] = 1.;
    a1[2] = 2.;
    a1[3] = 3.;
    hexed::Array<double> a2 = a0 + a1;
    REQUIRE(a0[1] == Catch::Approx(.1));
    REQUIRE(a1[1] == Catch::Approx(1.));
    REQUIRE(a2[1] == Catch::Approx(1.1));
    REQUIRE(a2[3] == Catch::Approx(3.3));
    REQUIRE_THROWS(a0(0) + a1);
    REQUIRE(a2.copy<int>()[2] == 2);
    hexed::Array<double> a3 = -a2;
    REQUIRE(a3[3] == Catch::Approx(-3.3));
    hexed::Array<double> a4 = 2.*a0 + 10.;
    REQUIRE(a4[1] == Catch::Approx(10.2));
    a4 = 3./a0/2.;
    REQUIRE(a4[2] == Catch::Approx(7.5));
  }
}
