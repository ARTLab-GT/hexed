#include <catch2/catch.hpp>

#include <Threaded_list.hpp>

TEST_CASE("Threaded_list")
{
  cartdg::Threaded_list<int> tl4 {4};
  REQUIRE(tl4.thread_begins().size() == 4);
  REQUIRE(tl4.thread_ends().size() == 4);
  SECTION("Construction edge cases")
  {
    cartdg::Threaded_list<int> tl1 {1};
    REQUIRE(tl1.thread_begins().size() == 1);
    REQUIRE(tl1.thread_ends().size() == 1);
    cartdg::Threaded_list<int> tln;
    unsigned n_begins = tln.thread_begins().size();
    REQUIRE(n_begins > 0);
    REQUIRE(tln.thread_ends().size() == n_begins);
    std::vector<cartdg::Threaded_list<int>> tls;
    REQUIRE_THROWS(tls.emplace_back(0));
  }
}
