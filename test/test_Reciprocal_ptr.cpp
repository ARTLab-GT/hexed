#include <string>
#include <catch2/catch_all.hpp>
#include <hexed/Reciprocal_ptr.hpp>

class Derived0 : public hexed::Mortal
{
  public:
  std::string s;
  Derived0(std::string s_) : s{s_} {}
};

class Derived1 : public hexed::Mortal
{
  public:
  int i;
  Derived1(int i_) : i{i_} {}
};

TEST_CASE("Reciprocal_ptr")
{
  Derived0 d00("ubiquitous");
  Derived1 d10(1903);
  hexed::Reciprocal_ptr<Derived0, Derived1> ptr0(&d00);
  hexed::Reciprocal_ptr<Derived1, Derived0> ptr1(&d10);

  REQUIRE(!ptr0);
  REQUIRE(!ptr0.paired());
  REQUIRE(ptr0.get() == nullptr);
  REQUIRE_THROWS(ptr0.value());

  ptr1.pair(ptr0);
  REQUIRE(ptr1);
  REQUIRE(ptr0);
  REQUIRE(ptr1.paired());
  REQUIRE(ptr0.paired());
  REQUIRE(ptr0.get() == &d10);
  REQUIRE(ptr1.get() == &d00);
  REQUIRE((*ptr0).i == 1903);
  REQUIRE(ptr1->s == "ubiquitous");
  REQUIRE(ptr1.value().s == "ubiquitous");

  ptr0.mine.set();
  REQUIRE(ptr0);
  REQUIRE(!ptr1);
  REQUIRE(ptr0.paired());
  REQUIRE(ptr1.paired());
  REQUIRE(ptr1.get() == nullptr);
  Derived0 d01("mendacious");
  ptr0.mine.set(&d01);
  REQUIRE(ptr1.get() == &d01);

  ptr0.unpair();
  REQUIRE(!ptr0);
  REQUIRE(!ptr1);
  REQUIRE(!ptr0.paired());
  REQUIRE(!ptr1.paired());
  ptr0.unpair();
  ptr0.pair(ptr1);

  {
    Derived0 d02("polyglottal");
    hexed::Reciprocal_ptr<Derived0, Derived1> ptr2(&d02);
    REQUIRE(ptr1.get() == &d01);
    ptr2.pair(ptr1);
    REQUIRE(!ptr0);
    REQUIRE(ptr1);
    REQUIRE(ptr2);
    REQUIRE(!ptr0.paired());
    REQUIRE(ptr1.paired());
    REQUIRE(ptr2.paired());
    REQUIRE(ptr0.get() == nullptr);
    REQUIRE(ptr1.get() == &d02);
    REQUIRE(ptr2.get() == &d10);
  }
  REQUIRE(!ptr1);
  REQUIRE(!ptr1.paired());
  REQUIRE(ptr1.get() == nullptr);

  ptr0.pair(ptr1);
  hexed::Reciprocal_ptr<Derived0, Derived1> ptr3(std::move(ptr0));
  REQUIRE(ptr3.mine.get() == &d01);
  REQUIRE(!ptr0);
  REQUIRE(ptr1);
  REQUIRE(ptr3);
  REQUIRE(!ptr0.paired());
  REQUIRE(ptr1.paired());
  REQUIRE(ptr3.paired());
  REQUIRE(ptr1.get() == &d01);
  REQUIRE(ptr3.get() == &d10);
}
