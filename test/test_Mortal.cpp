#include <catch2/catch_all.hpp>
#include <hexed/Mortal.hpp>

class Derived : public hexed::Mortal
{
  public:
  int i;
  Derived(int ii) : i{ii} {}
};

TEST_CASE("Mortal")
{
  Derived d0{1903};
  REQUIRE(d0.partners().size() == 0);
  hexed::Mortal_ptr<Derived> p0;
  REQUIRE(!p0);
  REQUIRE(p0.get() == nullptr);
  REQUIRE_THROWS(p0.value());
  hexed::Mortal_ptr<Derived> p1(&d0);
  REQUIRE(d0.partners().size() == 1);
  REQUIRE(p1);
  REQUIRE(p1.get() == &d0);
  REQUIRE(&(*p1) == &d0);
  REQUIRE(p1->i == 1903);
  REQUIRE(&p1.value() == &d0);
  p0.set(&d0);
  REQUIRE(d0.partners().size() == 2);
  REQUIRE(p0.get() == &d0);
  REQUIRE(p1.get() == &d0);

  {
    Derived d1{71};
    p1.set(&d1);
    REQUIRE(d0.partners().size() == 1);
    REQUIRE(p1.get() == &d1);
    REQUIRE(p0.get() == &d0);
  }
  REQUIRE(!p1);
  REQUIRE(p1.get() == nullptr);

  {
    hexed::Mortal_ptr<Derived> p2(&d0);
    REQUIRE(d0.partners().size() == 2);
  }
  REQUIRE(d0.partners().size() == 1);

  SECTION("Mortal(Mortal&&)") {
    REQUIRE(p0.get() == &d0);
    Derived d2 = std::move(d0);
    REQUIRE(d2.i == 1903);
    REQUIRE(d2.partners().size() == 1);
    REQUIRE(p0.get() == &d2);
    hexed::Mortal m;
    m = std::move(d2);
    REQUIRE_THROWS(p0.get());
  }

  SECTION("Mortal_ptr(Mortal_ptr&&)") {
    hexed::Mortal_ptr<Derived> p2 = std::move(p0);
    REQUIRE(d0.partners().size() == 1);
    REQUIRE(p2.get() == &d0);
    REQUIRE(!p0);
    REQUIRE(p0.get() == nullptr);
  }
}
