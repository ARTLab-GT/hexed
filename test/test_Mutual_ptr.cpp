#include <catch2/catch_all.hpp>
#include <hexed/Mutual_ptr.hpp>

TEST_CASE("Mutual_ptr")
{
  REQUIRE_THROWS(hexed::Mutual_ptr<double, int>(nullptr));
  double d = 8.314;
  int i = 6;
  hexed::Mutual_ptr<double, int> pdi(&d);
  hexed::Mutual_ptr<int, double> pid(&i);
  REQUIRE(pdi.mine() == Catch::Approx(8.314));
  REQUIRE(pid.mine() == 6);
  REQUIRE(!pdi);
  REQUIRE(!pid);
  REQUIRE(!pdi.partner());
  REQUIRE(!pdi.get());
  REQUIRE_THROWS(pdi.value());
  pdi.unpair();
  REQUIRE(!pdi);
  pid.pair(pdi);
  REQUIRE(pid);
  REQUIRE(pid.value() == Catch::Approx(8.314));
  REQUIRE(pdi.value() == 6);
  REQUIRE(*pdi == 6);
  REQUIRE(pdi.get() == &i);
  int* p = pdi;
  REQUIRE(p == &i);

  SECTION("unpairing")
  {
    pid.unpair();
    REQUIRE(!pid);
    REQUIRE(!pdi);
    REQUIRE(!pdi.get());
    REQUIRE_THROWS(pdi.value());
  }

  SECTION("pairing with an already paired ptr")
  {
    double e = 287.0528;
    hexed::Mutual_ptr<double, int> qdi(&e);
    qdi.pair(pid);
    REQUIRE(!pdi.partner());
    REQUIRE_THROWS(pdi.value());
    REQUIRE(!pdi);
    REQUIRE(qdi.value() == 6);
    REQUIRE(pid.value() == Catch::Approx(287.0528));
  }

  SECTION("move semantics")
  {
    hexed::Mutual_ptr<double, int> qdi(std::move(pdi));
    REQUIRE(qdi);
    REQUIRE(qdi.mine() == Catch::Approx(8.314));
    REQUIRE(qdi.value() == 6);
    qdi.mine() = 9.81;
    REQUIRE(pid.value() == Catch::Approx(9.81));
  }

  SECTION("destruction")
  {
    {
      int j = 7;
      hexed::Mutual_ptr<int, double> qid(&j);
      pdi.pair(qid);
      REQUIRE(!pid);
      REQUIRE(pdi.value() == 7);
      REQUIRE(qid);
    }
    REQUIRE(!pdi);
  }
}
