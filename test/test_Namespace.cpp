#include <catch2/catch_all.hpp>
#include <hexed/Namespace.hpp>

TEST_CASE("Namespace")
{
  auto space = std::make_shared<hexed::Namespace>();
  // assign some variables
  space->assign<int>("number0", 4);
  space->assign<double>("number1", 0.6);
  space->assign<std::string>("person", "Montgomery Knight");
  // check variable lookup mechanics
  REQUIRE(space->lookup<std::string>("person").value() == "Montgomery Knight");
  REQUIRE(!space->lookup<std::string>("unperson")); // "unperson" not assigned
  REQUIRE(!space->lookup<std::string>("number0")); // "number0" is not a string
  REQUIRE(space->lookup<double>("number1").value() == Catch::Approx(0.6));
  REQUIRE(space->lookup<double>("number0").value() == Catch::Approx(4.0)); // can read a double as an int
  REQUIRE(space->lookup<int>("number0").value() == 4);
  REQUIRE(!space->lookup<std::string>("number1")); // can't read a double as a string
  REQUIRE(space->exists("person"));
  REQUIRE(space->exists_recursive("person"));
  REQUIRE(!space->exists("unperson"));
  REQUIRE(!space->exists_recursive("unperson"));
  REQUIRE_THAT(space->names(), Catch::Matchers::UnorderedEquals(std::vector<std::string>{"number0", "number1", "person"}));
  // reassignment mechanics
  space->assign<std::string>("person", "Robert Jones");
  REQUIRE(space->lookup<std::string>("person").value() == "Robert Jones");
  REQUIRE_THROWS(space->assign<int>("person", 7)); // can't assign an int to a string
  space->assign<double>("number1", 3.14);
  REQUIRE(space->lookup<double>("number1").value() == Catch::Approx(3.14));
  space->assign<int>("number1", 7); // can assign an int to a double...
  REQUIRE(space->lookup<double>("number1").value() == Catch::Approx(7.0));
  REQUIRE(!space->lookup<int>("number1")); // ...but it's still a double
  space->assign_default<int>("number0", 6); // `assign_default` won't modify an existing value
  REQUIRE(space->lookup<int>("number0").value() == 4);
  space->assign_default<int>("number2", 6); // ...but it can create a new value
  REQUIRE(space->lookup<int>("number2").value() == 6);
  // read-only variables
  int _counter = 0;
  space->create("call_counter", new hexed::Namespace::Heisenberg<int>([&_counter]()->int{return _counter++;}));
  REQUIRE(space->lookup<int>("call_counter").value() == 0);
  REQUIRE(space->lookup<int>("call_counter").value() == 1);
  // super/sub namespaces
  auto sub = std::make_shared<hexed::Namespace>();
  sub->supers.push_back(space);
  auto subsub = std::make_shared<hexed::Namespace>();
  subsub->supers.push_back(sub);
  auto subsubsub = std::make_shared<hexed::Namespace>();
  subsubsub->supers.push_back(subsub);
  REQUIRE(sub->names().empty());
  REQUIRE(sub->lookup<int>("number0").value() == 4);
  REQUIRE(subsubsub->lookup<int>("number0").value() == 4);
  sub->assign<int>("number0", 6);
  REQUIRE(sub->lookup<int>("number0").value() == 6);
  REQUIRE(space->lookup<int>("number0").value() == 4);
  REQUIRE(!sub->exists("person"));
  REQUIRE(sub->exists_recursive("person"));
  sub->assign<int>("person", 0);
  REQUIRE(sub->lookup<int>("person").value() == 0);
  REQUIRE(!sub->lookup<std::string>("person"));
  REQUIRE(!sub->exists("unperson"));
  REQUIRE(!sub->exists_recursive("unperson"));
  // test basic array creation and lookup
  space->assign("data", hexed::Array<double>::make_uniform({6}, 287.0528));
  auto data0 = space->lookup<hexed::Array<double>>("data");
  REQUIRE(data0.value()[3] == Catch::Approx(287.0528));
  // test that lookups get you a reference rather than a copy
  data0.value()[1] = -1.;
  auto data1 = space->lookup<hexed::Array<double>>("data");
  REQUIRE(data1.value()[1] == Catch::Approx(-1.));
}
