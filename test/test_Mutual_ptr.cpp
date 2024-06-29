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

  SECTION("redundant pairing")
  {
    pdi.pair(*pdi.partner());
    REQUIRE(pdi);
    REQUIRE(pid);
    REQUIRE(pid.value() == Catch::Approx(8.314));
    REQUIRE(pdi.value() == 6);
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

#define REQ_CONTENTS(vec, ...) REQUIRE_THAT(vec, Catch::Matchers::UnorderedRangeEquals(__VA_ARGS__));

TEST_CASE("Multiple_ptr") {
  REQUIRE_THROWS(hexed::Multiple_ptr<int, std::string>(nullptr));
  int i = 43;
  std::string s0 = "ubiquitous";
  std::string s1 = "mendacious";
  std::string s2 = "polyglottal";
  std::unique_ptr<hexed::Multiple_ptr<int, std::string>> multi0(new hexed::Multiple_ptr<int, std::string>(&i));
  std::unique_ptr<hexed::Mutual_ptr<std::string, int>> mutual0(new hexed::Mutual_ptr<std::string, int>(&s0));
  std::unique_ptr<hexed::Mutual_ptr<std::string, int>> mutual1(new hexed::Mutual_ptr<std::string, int>(&s1));
  std::unique_ptr<hexed::Multiple_ptr<std::string, int>> multi1(new hexed::Multiple_ptr<std::string, int>(&s2));

  SECTION("connecting and disconnecting") {
    multi0->add(*mutual0);
    mutual1->pair(*multi0);
    multi1->add(*multi0);
    REQUIRE(&multi0->mine() == &i);
    REQ_CONTENTS(multi0->partners(), std::vector<void*>{mutual0.get(), mutual1.get(), multi1.get()});
    REQ_CONTENTS(multi1->partners(), std::vector<void*>{multi0.get()});
    REQUIRE(mutual0->partner() == multi0.get());
    std::unique_ptr<hexed::Multiple_ptr<int, std::string>> multi2(new hexed::Multiple_ptr<int, std::string>(&i));
    multi2->add(*mutual0);
    REQ_CONTENTS(multi0->partners(), std::vector<void*>{mutual1.get(), multi1.get()});
    REQ_CONTENTS(multi2->partners(), std::vector<void*>{mutual0.get()});
    multi0->remove(*mutual1);
    REQ_CONTENTS(multi0->partners(), std::vector<void*>{multi1.get()});
    REQUIRE(!*mutual1);
    multi1.reset();
    multi2.reset();
    REQUIRE(!*mutual0);
    REQUIRE(multi0->partners().size() == 0);
  }

  SECTION("move semantics") {
    int j = 0;
    hexed::Multiple_ptr<int, std::string> multi3(&i);
    hexed::Multiple_ptr<int, std::string> multi4(&j);
    hexed::Multiple_ptr<std::string, int> multi5(&s1);
    mutual0->pair(multi3);
    multi3.add(*multi1);
    multi4.add(multi5);
    multi4 = std::move(multi3);
    REQUIRE(mutual0->partner() == &multi4);
    REQ_CONTENTS(multi4.partners(), std::vector<void*>{mutual0.get(), multi1.get()});
    REQUIRE(&multi4.mine() == &i);
    REQUIRE(multi3.partners().empty());
    hexed::Multiple_ptr<int, std::string> multi6(std::move(multi4));
    REQUIRE(&multi4.mine() == &i);
    REQ_CONTENTS(multi6.partners(), std::vector<void*>{mutual0.get(), multi1.get()});
    REQUIRE(multi4.partners().empty());
  }
}
