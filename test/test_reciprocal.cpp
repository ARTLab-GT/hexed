#include <string>
#include <catch2/catch_all.hpp>
#include <hexed/reciprocal.hpp>

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

TEST_CASE("reciprocal")
{
  Derived0 d00("ubiquitous");
  Derived0 d01("mendacious");
  Derived1 d10(1903);
  hexed::Reciprocal_ptr<Derived0, Derived1> ptr0(&d00);
  hexed::Reciprocal_ptr<Derived1, Derived0> ptr1(&d10);

  SECTION("Reciprocal_ptr") {
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

  SECTION("Reciprocal_list") {
    hexed::Reciprocal_list<Derived0, Derived1> list0(&d01);
    hexed::Reciprocal_list<Derived1, Derived0> list1(&d10);
    #define REQ_SAME_ADDRS(refs, T, U, ...) { \
        auto addrs_seq = refs.transform<void*>([](hexed::mutual::Base<T, U>& ref)->void*{return &ref;}); \
        std::vector<void*> addrs(addrs_seq.begin(), addrs_seq.end()); \
        REQUIRE_THAT(addrs, Catch::Matchers::UnorderedRangeEquals(std::vector<void*>__VA_ARGS__)); \
      }
    REQUIRE(!ptr0);
    REQUIRE(!list0.partners());
    REQUIRE(!list1.partners());
    list1.add(list0);
    ptr0.pair(list1);
    REQ_SAME_ADDRS(list0.partners(), Derived1, Derived0, {&list1});
    REQ_SAME_ADDRS(list1.partners(), Derived0, Derived1, {&ptr0, &list0});
    REQUIRE_THAT(list1.theirs(), Catch::Matchers::UnorderedRangeEquals(std::vector<void*>{&d00, &d01}));
    REQUIRE(ptr0.get() == &d10);

    SECTION("unpair/remove") {
      list0.remove(ptr1);
      REQ_SAME_ADDRS(list0.partners(), Derived1, Derived0, {&list1});
      ptr0.unpair();
      list0.remove(list1);
      REQUIRE(list0.partners().empty());
      REQ_SAME_ADDRS(list1.partners(), Derived0, Derived1, {&ptr0});
    }

    SECTION("clear") {
      list1.clear();
      REQUIRE(!ptr0);
      REQUIRE(list1.partners().empty());
      REQUIRE(list0.partners().empty());
    }

    SECTION("multiple add") {
      list0.add(list1);
      REQ_SAME_ADDRS(list0.partners(), Derived1, Derived0, {&list1});
      REQ_SAME_ADDRS(list1.partners(), Derived0, Derived1, {&ptr0, &list0});
    }
  }
}
