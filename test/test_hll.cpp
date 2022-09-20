#include <catch2/catch.hpp>
#include <hexed/hll.hpp>

TEST_CASE("hll")
{
  double mass = 1.225;
  double pressure = 101325;
  double read[20] {};
  double backup[20] {};

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
    hexed::hll<3, 2>(read, 0, 1.4);
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
    hexed::hll<3, 2>(read, 1, 1.4);
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
    hexed::hll<3, 2>(read, 0, 1.4);
    REQUIRE(read[2*3     ] == Approx(0).margin(1e-8));
    REQUIRE(read[2*3 + 10] == Approx(0).margin(1e-8));
  }
}
