#include <catch2/catch_all.hpp>
#include <hexed/hdf5_utils.hpp>

TEST_CASE("hdf5_utils") {
  {
    H5::H5File file("hdf5_test.h5", H5F_ACC_TRUNC);
    hexed::hdf5_utils::add_attr(file, "ubiquitous", 6);
    hexed::hdf5_utils::add_attr<hexed::Int>(file, "mendacious", 7);
    hexed::hdf5_utils::add_attr(file, "polyglottal", 0.42);
    REQUIRE(hexed::hdf5_utils::get_attr<int>(file, "ubiquitous") == 6);
    REQUIRE(hexed::hdf5_utils::get_attr<hexed::Int>(file, "mendacious") == 7);
    REQUIRE(hexed::hdf5_utils::get_attr<double>(file, "polyglottal") == Catch::Approx(0.42));
  }
  REQUIRE(hexed::hdf5_utils::get_attr<int>("hdf5_test.h5", "ubiquitous") == 6);
  REQUIRE(hexed::hdf5_utils::get_attr<hexed::Int>("hdf5_test.h5", "mendacious") == 7);
  REQUIRE(hexed::hdf5_utils::get_attr<double>("hdf5_test.h5", "polyglottal") == Catch::Approx(0.42));
}
