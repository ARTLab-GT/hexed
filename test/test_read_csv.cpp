#include <catch2/catch_all.hpp>
#include <hexed/read_csv.hpp>

TEST_CASE("read_csv")
{
  auto data = hexed::read_csv(hexed::config::root_dir + std::string("/test/test_csv.csv"));
  REQUIRE(data.rows() == 2);
  REQUIRE(data.cols() == 3);
  Eigen::MatrixXd correct(2, 3);
  correct << 0.1, 1e3, 10., -.3, 4., 1.;
  REQUIRE((data - correct).norm() < 1e-12);
}
