#include <thread>
#include <fstream>
#include <catch2/catch.hpp>
#include <hexed/Stopwatch_tree.hpp>

TEST_CASE("Stopwatch_tree")
{
  // no good way to really test this, so just do a silly example and look at the output
  hexed::Stopwatch_tree rocinante {"gunship"};
  rocinante.children.emplace("armaments", "weapon system");
  rocinante.children.at("armaments").children.emplace("defensive", "PDC");
  rocinante.children.at("armaments").children.emplace("long range", "missile");
  rocinante.children.at("armaments").children.emplace("close range", "railgun");
  rocinante.children.emplace("fuel", "kg");
  rocinante.children.emplace("equipment", "EVA suit");
  rocinante.children.at("equipment").children.emplace("power armor", "set");

  rocinante.stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  rocinante.children.at("armaments").stopwatch.start();

  rocinante.children.at("armaments").children.at("defensive").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  rocinante.children.at("armaments").children.at("defensive").stopwatch.pause();
  rocinante.children.at("armaments").children.at("defensive").work_units_completed += 5;

  rocinante.children.at("armaments").children.at("long range").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  rocinante.children.at("armaments").children.at("long range").stopwatch.pause();
  rocinante.children.at("armaments").children.at("long range").work_units_completed += 10;

  rocinante.children.at("armaments").children.at("close range").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(20));
  rocinante.children.at("armaments").children.at("close range").stopwatch.pause();
  rocinante.children.at("armaments").children.at("close range").work_units_completed += 1;

  rocinante.children.at("armaments").stopwatch.pause();
  rocinante.children.at("armaments").work_units_completed += 3;

  rocinante.children.at("fuel").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  rocinante.children.at("fuel").stopwatch.pause();
  rocinante.children.at("fuel").work_units_completed += 1000;

  rocinante.children.at("equipment").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
  rocinante.children.at("equipment").stopwatch.pause();
  rocinante.children.at("equipment").work_units_completed += 5;
  rocinante.children.at("equipment").children.at("power armor").stopwatch.start();
  std::this_thread::sleep_for(std::chrono::milliseconds(0));
  rocinante.children.at("equipment").children.at("power armor").stopwatch.pause();
  rocinante.children.at("equipment").children.at("power armor").work_units_completed += 1;

  rocinante.stopwatch.pause();
  rocinante.work_units_completed += 1;
  std::ofstream output ("build_a_spaceship.txt");
  output << rocinante.report();
}
