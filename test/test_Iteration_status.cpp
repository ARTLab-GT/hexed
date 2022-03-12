#include <catch2/catch.hpp>
#include <Iteration_status.hpp>

bool in(std::string str0, std::string str1)
{
  return str1.find(str0) != std::string::npos;
}

TEST_CASE("Iteration_status")
{
  cartdg::Iteration_status stat;
  std::vector<std::string> lines;
  lines.push_back(stat.header());
  stat.art_visc_iters = 3;
  lines.push_back(stat.report());
  stat.flow_time = 2.314;
  stat.time_step = 0.7;
  stat.iteration = int(1e6);
  lines.push_back(stat.report());
  stat.flow_time = 2.3141;
  stat.time_step = -0.7;
  ++stat.iteration;
  lines.push_back(stat.report());
  // test, as well as possible, that the strings contain the information they are supposed to
  REQUIRE(in("iteration", lines[0]));
  REQUIRE(in("flow time", lines[0]));
  REQUIRE(in("time step", lines[0]));
  REQUIRE(in("3", lines[1]));
  REQUIRE(in("1000000", lines[2]));
  REQUIRE(in("2", lines[2]));
  REQUIRE(in("3141", lines[3]));
  REQUIRE(in(" -", lines[3]));
  REQUIRE(in("7", lines[3]));
  // the values should be neatly aligned, so the strings should all be the same length
  REQUIRE(lines[0].size() == lines[1].size());
  REQUIRE(lines[0].size() == lines[2].size());
  // trailing spaces are annoying
  for (auto line : lines) REQUIRE(line.back() != ' ');
}
