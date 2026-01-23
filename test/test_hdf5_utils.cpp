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
    hsize_t dims[2];
    dims[0] = 3;
    dims[1] = 2;
    auto dset = file.createDataSet("/some_data", hexed::hdf5_utils::type<double>(), H5::DataSpace(2, dims));
    hexed::hdf5_utils::write(dset, 0, 0, 0.0);
    hexed::hdf5_utils::write(dset, 0, 1, 0.1);
    hexed::hdf5_utils::write(dset, 1, 0, 0.2);
    hexed::hdf5_utils::write(dset, 1, 1, 0.3);
    hexed::hdf5_utils::write(dset, 2, 0, 0.4);
    hexed::hdf5_utils::write(dset, 2, 1, 0.5);
    REQUIRE(hexed::hdf5_utils::read<double>(dset, 0, 1) == Catch::Approx(.1));
    REQUIRE(hexed::hdf5_utils::read<double>(dset, 1, 0) == Catch::Approx(.2));
    REQUIRE(hexed::hdf5_utils::read<double>(dset, 2, 1) == Catch::Approx(.5));
  }
  REQUIRE(hexed::hdf5_utils::get_attr<int>("hdf5_test.h5", "ubiquitous") == 6);
  REQUIRE(hexed::hdf5_utils::get_attr<hexed::Int>("hdf5_test.h5", "mendacious") == 7);
  REQUIRE(hexed::hdf5_utils::get_attr<double>("hdf5_test.h5", "polyglottal") == Catch::Approx(0.42));
}
