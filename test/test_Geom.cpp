#include <catch2/catch_all.hpp>
#include <hexed/cad_geom.hpp>

TEST_CASE("Geom") {
  hexed::cad_geom::Geom geom("../test_assets/cylinder_extruded.iges");
  geom.visualize("cylinder_extruded");
}
