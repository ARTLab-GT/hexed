#include <catch2/catch.hpp>
#include <hexed/erase_if.hpp>

bool is_6(const int& i)
{
  return i == 6;
}

TEST_CASE("erase_if")
{
  std::vector<int> small {1, 6, 3, 4, 6, 6, 0, 6};
  std::vector<int> big (100, 7);
  big[0] = 6;
  big[4] = 6;
  big[99] = 6;
  hexed::erase_if(small, &is_6);
  hexed::erase_if(big, &is_6);
  std::vector<int> correct_small {1, 3, 4, 0};
  REQUIRE(std::equal(correct_small.begin(), correct_small.end(), small.begin()));
  REQUIRE(big.size() == 97);
}
