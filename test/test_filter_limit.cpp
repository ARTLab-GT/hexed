#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/math.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/filter_limit.hpp>

TEST_CASE("filter_limit")
{
  int row_size = hexed::config::max_row_size;
  Eigen::VectorXd data(row_size*row_size);
  hexed::Gauss_legendre basis(row_size);
  // set polynomial to interpolation of smooth transcendental function
  for (int i = 0; i < row_size; ++i) {
    for (int j = 0; j < row_size; ++j) {
      data(i*row_size + j) = std::exp(basis.node(i) + 0.5*basis.node(j));
    }
  }
  auto norm = [&]() {return data.dot(hexed::math::pow_outer(basis.node_weights(), 2).cwiseProduct(data));};
  double before = norm();
  double decay_rate = 0.2;
  hexed::filter_limit(2, data.data(), basis, decay_rate);
  // check that the function has not been modified much
  CHECK(norm()/before > 0.9);
  srand(406);
  // set to a very not-smooth function
  for (int i = 0; i < row_size; ++i) {
    for (int j = 0; j < row_size; ++j) {
      data(i*row_size + j) = 2e-3*(rand()%1000) - 1;
    }
  }
  before = norm();
  hexed::filter_limit(2, data.data(), basis, decay_rate);
  // check that the function has been modified significantly
  CHECK(norm()/before < 0.1);
}
