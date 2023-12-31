#include <catch2/catch_all.hpp>
#include <Eigen/Dense>
#include <hexed/Physical_basis.hpp>

TEST_CASE("Physical_basis")
{
  srand(3);

  SECTION("input size checking")
  {
    std::vector<double> points (50);
    std::vector<hexed::Physical_basis> bases;
    REQUIRE_THROWS(bases.emplace_back(1, 5, points));
    REQUIRE_THROWS(bases.emplace_back(2, 4, points));
    REQUIRE_THROWS(bases.emplace_back(2, 6, points));
    bases.emplace_back(2, 5, points);
    std::vector<double> short_points (4*16);
    REQUIRE_THROWS(bases.emplace_back(4, 2, short_points));
  }

  SECTION("basis size calculation")
  {
    {
      std::vector<double> points(2);
      hexed::Physical_basis basis {2, 1, points};
      REQUIRE(basis.size() == 1);
    }
    {
      std::vector<double> points(8);
      hexed::Physical_basis basis {2, 2, points};
      REQUIRE(basis.size() == 3);
    }
    {
      std::vector<double> points(32);
      hexed::Physical_basis basis {2, 4, points};
      REQUIRE(basis.size() == 10);
    }
    {
      std::vector<double> points(3);
      hexed::Physical_basis basis {3, 1, points};
      REQUIRE(basis.size() == 1);
    }
    {
      std::vector<double> points(24);
      hexed::Physical_basis basis {3, 2, points};
      REQUIRE(basis.size() == 4);
    }
    {
      std::vector<double> points(192);
      hexed::Physical_basis basis {3, 4, points};
      REQUIRE(basis.size() == 20);
    }
  }

  SECTION("numerics")
  {
    int n_qpoint = 27;
    std::vector<double> points (3*n_qpoint);
    Eigen::MatrixXd values (n_qpoint, 2);
    Eigen::VectorXd weights (n_qpoint);
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        // random points within a 0.1 by 0.1 by 0.1 box at 100., 100., 100.
        points[i_dim*n_qpoint + i_qpoint] = (rand() % 1000)/10000. + 100.;
      }
      for (int i_poly = 0; i_poly < 2; ++i_poly) {
        values(i_qpoint, i_poly) = (rand() % 1000)/500. - 1.; // polynomial value between -1. and 1.
      }
      weights(i_qpoint) = (rand() % 1000)/2000. + 0.5; // weights between 0.5 and 1.
    }
    hexed::Physical_basis basis (3, 3, points);
    Eigen::MatrixXd mat (n_qpoint, basis.size());
    for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
      for (int i_basis = 0; i_basis < basis.size(); ++i_basis) {
        mat(i_qpoint, i_basis) = basis.evaluate(i_qpoint, i_basis);
      }
    }

    SECTION("evaluation")
    {
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr (mat);
      REQUIRE(qr.isInjective()); // check that the basis polynomials are linearly independent
      Eigen::VectorXd values;
      // check that the constant function 1 is in the space
      values = Eigen::VectorXd::Ones(n_qpoint);
      REQUIRE((mat*qr.solve(values) - values).norm() < 1e-8);
      // check that ((x_i - 100)/0.1)^2 is in the space
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
          double coord = (points[i_dim*n_qpoint + i_qpoint] - 100.)/0.1;
          values[i_qpoint] = coord*coord;
        }
        REQUIRE((mat*qr.solve(values) - values).norm() < 1e-8);
      }
      // check that ((x_0 - 100)/0.1)*((x_1 - 100)/0.1) is in the space
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        values[i_qpoint] = 1.;
        for (int i_dim = 0; i_dim < 2; ++i_dim) {
          double coord = (points[i_dim*n_qpoint + i_qpoint] - 100.)/0.1;
          values[i_qpoint] *= coord;
        }
      }
      REQUIRE((mat*qr.solve(values) - values).norm() < 1e-8);
      // check that ((x_0 - 100)/0.1)*((x_0 - 100)/0.1)*((x_2 - 100)/0.1) is NOT in the space
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        values[i_qpoint] = 1.;
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
          double coord = (points[i_dim*n_qpoint + i_qpoint] - 100.)/0.1;
          values[i_qpoint] *= coord;
        }
      }
      REQUIRE((mat*qr.solve(values) - values).norm() > 1e-3);
    }

    SECTION("projection")
    {
      Eigen::MatrixXd projected = basis.projection(values, weights);
      Eigen::MatrixXd twice = basis.projection(projected, weights);
      REQUIRE((projected - twice).norm() < 1e-8);
      REQUIRE((mat.transpose()*weights.asDiagonal()*(projected - values)).norm() < 1e-8);
    }
  }
}
