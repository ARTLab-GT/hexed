#include <catch2/catch_all.hpp>
#include <hexed/Neighbor_connection.hpp>

TEST_CASE("Neighbor_connection") {
  hexed::Storage_params params {2, 4, 2, 2};
  hexed::Face f0(params, 0, 0, false);
  hexed::Face f1(params, 1, 0, true);
  std::unique_ptr<hexed::Face> f2 {new hexed::Face(params, 1, 1, true)};
  {
    hexed::Neighbor_connection con(params, {&f0, f2.get()});
    REQUIRE(con.alive());
    REQUIRE(&con.face(0) == &f0);
    REQUIRE(&con.face(1) == f2.get());
    REQUIRE(f0.neighbor_connection() == &con);
    REQUIRE(f2->neighbor_connection() == &con);
    REQUIRE(f0.connected());
    REQUIRE(f2->connected());
    REQUIRE(&con.opposite_face(*f2) == &f0);
    REQUIRE(&con.opposite_face(f0) == f2.get());
    REQUIRE(con.is_deformed() == false);
    f2.reset();
    REQUIRE(!con.alive());
    REQUIRE(&con.face(0) == &f0);
    REQUIRE_THROWS(con.face(1));
    REQUIRE_THROWS(hexed::Neighbor_connection (params, {&f0, &f1}));
    REQUIRE(con.is_deformed() == false);
  }
  hexed::Neighbor_connection con(params, {&f0, &f1});
  REQUIRE(con.alive());
  f0.disconnect();
  REQUIRE(!f0.connected());
  REQUIRE(!con.alive());
}
