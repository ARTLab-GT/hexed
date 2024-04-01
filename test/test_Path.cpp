#include <catch2/catch_all.hpp>
#include <hexed/Path.hpp>
#include <fstream>

TEST_CASE("Path")
{
  hexed::Path p("lib");
  REQUIRE(!p.find("libhdf5_cpp.so").empty()); // this should exist somewhere, since HDF5 is a mandatory dependency
  // find a file in the home directory
  std::ofstream file("find_this");
  file.close();
  REQUIRE(!p.find("find_this").empty());
  REQUIRE(p.find("hexed_this_does_not_exist").empty());
  REQUIRE(!p.find("/dev/null").empty());
}
