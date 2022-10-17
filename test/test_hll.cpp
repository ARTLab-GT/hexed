#include <catch2/catch.hpp>
#include <Eigen/Dense>
#include <hexed/hll.hpp>

TEST_CASE("hll::inviscid")
{
  double mass = 1.225;
  double pressure = 101325;
  double read[20] {};
  double backup[20] {};
  Eigen::MatrixXd normal;

  SECTION("Reasonable flow")
  {
    double velocity0 [] {3*340, 2*340};
    double velocity1 [] {-2*340, -3*340};
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        backup[i + 10*j + 0] = mass*velocity0[j];
        backup[i + 10*j + 2] = mass*velocity1[j];
        backup[i + 10*j + 4] = 0;
        backup[i + 10*j + 6] = mass;
        backup[i + 10*j + 8] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                        + velocity1[j]*velocity1[j]);
      }
    }
    for (int i = 0; i < 20; ++i) read[i] = backup[i];
    normal = Eigen::VectorXd::Ones(2)*Eigen::VectorXd::Unit(3, 0).transpose();
    hexed::hll::inviscid<3, 2>(read, normal.data(), 1.4);
    for (int j = 0; j < 2; ++j) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        double correct_flux = backup[2*i_var];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= 3*340;
        if (i_var == 0) correct_flux += pressure;
        REQUIRE(read[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(read[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) read[i] = backup[i];
    normal = Eigen::VectorXd::Ones(2)*Eigen::VectorXd::Unit(3, 1).transpose();
    hexed::hll::inviscid<3, 2>(read, normal.data(), 1.4);
    for (int j = 0; j < 2; ++j) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        double correct_flux = backup[2*i_var + 10];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= -3*340;
        if (i_var == 1) correct_flux += pressure;
        REQUIRE(read[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(read[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) read[i] = backup[i];
    Eigen::Vector3d n = Eigen::Vector3d{.8, .1, 0.};
    normal = Eigen::VectorXd::Ones(2)*n.transpose();
    hexed::hll::inviscid<3, 2>(read, normal.data(), 1.4);
    for (int j = 0; j < 2; ++j) {
      for (int i_var = 0; i_var < 5; ++i_var) {
        double correct_flux = backup[2*i_var];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= (3*n[0] - 2*n[1])*340;
        if (i_var < 3) correct_flux += pressure*n[i_var];
        REQUIRE(read[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(read[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
      }
    }
  }

  SECTION("Opposing supersonic flows")
  {
    double velocity0 [] {680, -680};
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        read[i + 10*j + 0] = mass*velocity0[j];
        read[i + 10*j + 2] = 0;
        read[i + 10*j + 4] = 0;
        read[i + 10*j + 6] = mass;
        read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*velocity0[j]*velocity0[j];
      }
    }
    normal = Eigen::VectorXd::Unit(3, 0)*Eigen::VectorXd::Ones(2).transpose();
    hexed::hll::inviscid<3, 2>(read, normal.data(), 1.4);
    REQUIRE(read[2*3     ] == Approx(0).margin(1e-8));
    REQUIRE(read[2*3 + 10] == Approx(0).margin(1e-8));
  }
}

TEST_CASE("hll::advection")
{
  double state [2][4][3];
  double nrml [2][3];
  double scalar [2] {1.3, 0.7};
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (int i_qpoint = 0; i_qpoint < 3; ++i_qpoint) {
      // set surface normal
      nrml[0][i_qpoint] = .4;
      nrml[1][i_qpoint] = .5;
      // set advection velocity
      state[i_side][0][i_qpoint] = -.1;
      state[i_side][1][i_qpoint] = -.9;
      // set scalar state;
      state[i_side][2][i_qpoint] = scalar[i_side];
    }
  }
  hexed::hll::advection<2, 3>(state[0][0], nrml[0]);
  for (int i_side = 0; i_side < 2; ++i_side) {
    for (int i_qpoint = 0; i_qpoint < 3; ++i_qpoint) {
      REQUIRE(state[i_side][0][i_qpoint] == Approx(-.1));
      REQUIRE(state[i_side][1][i_qpoint] == Approx(-.9));
      REQUIRE(state[i_side][2][i_qpoint] == Approx((-.1*.4 - .9*.5)*0.7));
    }
  }
}
