#include <catch.hpp>

#include <kernels/neighbor/read_copy.hpp>
#include <kernels/neighbor/write_copy.hpp>
#include <kernels/neighbor/cpg_euler_hll.hpp>
#include <kernels/neighbor/cpg_euler_hll_deformed.hpp>

TEST_CASE("neighbor kernel read_copy<>()")
{
  SECTION("1D")
  {
    double read [12];
    double write [3];
    for (int i = 0; i < 12; ++i) read[i] = i/10.;
    cartdg::read_copy<3, 4, 4>(&read[0], &write[0], 1, 0);
    REQUIRE(write[0] == .0);
    REQUIRE(write[1] == .4);
    REQUIRE(write[2] == .8);

    cartdg::read_copy<3, 4, 4>(&read[0], &write[0], 1, 1);
    REQUIRE(write[0] == .3);
    REQUIRE(write[1] == .7);
    REQUIRE(write[2] == 1.1);
  }

  SECTION("2D")
  {
    double read[18];
    double write[6];
    for (int i = 0; i < 18; ++i) read[i] = i/10.;
    cartdg::read_copy<2, 9, 3>(&read[0], &write[0], 3, 0);
    REQUIRE(write[0] == .0);
    REQUIRE(write[1] == .1);
    REQUIRE(write[2] == .2);
    REQUIRE(write[3] == .9);
    REQUIRE(write[4] == 1.0);
    REQUIRE(write[5] == 1.1);

    cartdg::read_copy<2, 9, 3>(&read[0], &write[0], 3, 1);
    REQUIRE(write[0] == .6);
    REQUIRE(write[1] == .7);
    REQUIRE(write[3] == 1.5);

    cartdg::read_copy<2, 9, 3>(&read[0], &write[0], 1, 1);
    REQUIRE(write[0] == .2);
    REQUIRE(write[1] == .5);
    REQUIRE(write[2] == .8);
    REQUIRE(write[5] == 1.7);
  }

  SECTION("3D")
  {
    double read[27];
    double write[9];
    for (int i = 0; i < 27; ++i) read[i] = i/10.;
    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 9, 0);
    REQUIRE(write[0] == .0);
    REQUIRE(write[1] == .1);
    REQUIRE(write[2] == .2);
    REQUIRE(write[3] == .3);
    REQUIRE(write[8] == .8);
    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 9, 1);
    REQUIRE(write[0] == 1.8);
    REQUIRE(write[7] == 2.5);
    REQUIRE(write[8] == 2.6);

    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 3, 0);
    REQUIRE(write[0] == .0);
    REQUIRE(write[1] == .1);
    REQUIRE(write[3] == .9);
    REQUIRE(write[8] == 2.);
    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 3, 1);
    REQUIRE(write[0] == .6);
    REQUIRE(write[2] == .8);
    REQUIRE(write[8] == 2.6);

    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 1, 0);
    REQUIRE(write[0] == .0);
    REQUIRE(write[1] == .3);
    REQUIRE(write[8] == 2.4);
    cartdg::read_copy<1, 27, 3>(&read[0], &write[0], 1, 1);
    REQUIRE(write[0] == .2);
    REQUIRE(write[7] == 2.3);
    REQUIRE(write[8] == 2.6);
  }
}

TEST_CASE("neighbor kernel cartdg::write_copy<>()")
{
  SECTION("1D")
  {
    double read [12];
    double write0 [3];
    double write1 [12];
    for (int i = 0; i < 12; ++i)
    {
      read[i] = i/10.;
      write1[i] = 0.;
    }
    #define COPY(stride, ipf) \
      cartdg::read_copy <3, 4, 4>(&read  [0], &write0[0], stride, ipf); \
      cartdg::write_copy<3, 4, 4>(&write0[0], &write1[0], stride, ipf);
    COPY(3, 0);
    COPY(3, 1);
    COPY(1, 0);
    COPY(1, 1);
    #undef COPY
    for (int i = 0; i < 12; ++i)
    {
      if ((i%4 == 0) || (i%4 == 3))
      {
        REQUIRE(read[i] == write1[i]);
      }
      else
      {
        REQUIRE(write1[i] == 0);
      }
    }
  }

  SECTION("2D")
  {
    double read[18];
    double write0[6];
    double write1[18];
    int n_copies[9] {2, 1, 2, 1, 0, 1, 2, 1, 2};
    for (int i = 0; i < 18; ++i)
    {
      read[i] = i/10.;
      write1[i] = 0.;
    }
    #define COPY(stride, ipf) \
      cartdg::read_copy <2, 9, 3>(&read  [0], &write0[0], stride, ipf); \
      cartdg::write_copy<2, 9, 3>(&write0[0], &write1[0], stride, ipf);
    COPY(3, 0);
    COPY(3, 1);
    COPY(1, 0);
    COPY(1, 1);
    #undef COPY
    for (int i = 0; i < 18; ++i)
    {
      REQUIRE(n_copies[i%9]*read[i] == write1[i]);
    }
  }

  SECTION("3D")
  {
    double read[27];
    double write0[9];
    double write1[27];
    int n_copies_2d[9] {2, 1, 2, 1, 0, 1, 2, 1, 2};
    int n_copies[27];
    for (int i = 0; i < 9; ++i)
    {
      n_copies[i + 9] = n_copies_2d[i];
      n_copies[i + 18] = n_copies[i] = n_copies_2d[i] + 1;
    }
    for (int i = 0; i < 27; ++i)
    {
      read[i] = i/10.;
      write1[i] = 0.;
    }
    #define COPY(stride, ipf) \
      cartdg::read_copy <1, 27, 3>(&read  [0], &write0[0], stride, ipf); \
      cartdg::write_copy<1, 27, 3>(&write0[0], &write1[0], stride, ipf);
    COPY(9, 0);
    COPY(9, 1);
    COPY(3, 0);
    COPY(3, 1);
    COPY(1, 0);
    COPY(1, 1);
    #undef COPY
    for (int i = 0; i < 27; ++i)
    {
      REQUIRE(n_copies[i]*read[i] == write1[i]);
    }
  }

}

TEST_CASE("cpg_euler_hll")
{
  double mass = 1.225;
  double pressure = 101325;
  double read[20] {};
  double write[20] {};
  double mult = 0.7;
  SECTION("Reasonable flow")
  {
    double velocity0 [] {3*340, 2*340};
    double velocity1 [] {-2*340, -3*340};
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        read[i + 10*j + 0] = mass*velocity0[j];
        read[i + 10*j + 2] = mass*velocity1[j];
        read[i + 10*j + 4] = 0;
        read[i + 10*j + 6] = mass;
        read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                      + velocity1[j]*velocity1[j]);
      }
    }
    cartdg::cpg_euler_hll<3, 2>(&read[0], &write[0], mult, 0, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_d_flux = 3*340*read[2*i_var] - 2*340*read[2*i_var + 10];
        if (i_var == 4) correct_d_flux += 340*pressure;
        correct_d_flux *= j;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(0.7*correct_d_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(0.7*correct_d_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) write[i] = 0;
    cartdg::cpg_euler_hll<3, 2>(&read[0], &write[0], mult, 1, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_d_flux = 3*340*read[2*i_var + 10] - 2*340*read[2*i_var];
        if (i_var == 4) correct_d_flux += 340*pressure;
        correct_d_flux *= (1 - j);
        REQUIRE(write[10*j + 2*i_var    ] == Approx(0.7*correct_d_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(0.7*correct_d_flux).margin(1.e-8));
      }
    }
  }

  SECTION("Opposing supersonic flows")
  {
    double velocity0 [] {680, -680};
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        read[i + 10*j + 0] = mass*velocity0[j];
        read[i + 10*j + 2] = 0;
        read[i + 10*j + 4] = 0;
        read[i + 10*j + 6] = mass;
        read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*velocity0[j]*velocity0[j];
      }
    }

    for (int i = 0; i < 20; ++i) write[i] = 0;
    cartdg::cpg_euler_hll<3, 2>(&read[0], &write[0], mult, 0, 1.4);
    REQUIRE(write[2*3     ] == Approx(0.7*mass*680));
    REQUIRE(write[2*3 + 10] == Approx(0.7*mass*680));
  }
}

TEST_CASE("cpg_euler_hll_deformed")
{
  double mass = 1.225;
  double pressure = 101325;
  double energy = pressure/0.4 + 0.5*mass*680*680;
  double read[20] {};
  double write[20] {};
  double mult = 0.7;
  bool flip [] {false, false};

  SECTION("Reasonable flow")
  {
    double jacobian [2][3][3][2] {};
    double velocity0 [] {3*340, 2*340};
    double velocity1 [] {-2*340, -3*340};
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        read[i + 10*j + 0] = mass*velocity0[j];
        read[i + 10*j + 2] = mass*velocity1[j];
        read[i + 10*j + 4] = 0;
        read[i + 10*j + 6] = mass;
        read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                      + velocity1[j]*velocity1[j]);
        for (int k = 0; k < 3; ++k)
        {
          jacobian[j][k][k][i] = 1.;
        }
      }
    }
    {
      int i_axis_arg [] {0, 0};
      cartdg::cpg_euler_hll_deformed<3, 2>(&read[0], &write[0],
                                           &(jacobian[0][0][0][0]), mult, i_axis_arg, flip, 1.4);
    }
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_d_flux = 3*340*read[2*i_var] - 2*340*read[2*i_var + 10];
        if (i_var == 4) correct_d_flux += 340*pressure;
        correct_d_flux *= j;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(0.7*correct_d_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(0.7*correct_d_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) write[i] = 0;
    {
      int i_axis_arg [] {1, 1};
      cartdg::cpg_euler_hll_deformed<3, 2>(&read[0], &write[0],
                                           &(jacobian[0][0][0][0]), mult, i_axis_arg, flip, 1.4);
    }
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_d_flux = 3*340*read[2*i_var + 10] - 2*340*read[2*i_var];
        if (i_var == 4) correct_d_flux += 340*pressure;
        correct_d_flux *= (1 - j);
        REQUIRE(write[10*j + 2*i_var    ] == Approx(0.7*correct_d_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(0.7*correct_d_flux).margin(1.e-8));
      }
    }
  }

  SECTION("normal flipping")
  {
    for (int i_give : {0, 1})
    {
      double jacobian [2][3][3][2] {};
      double velocity0 [2] {};
      double velocity1 [2] {};
      velocity1[i_give] = 3*340;
      velocity1[1 - i_give] = 2*340;
      for (int i = 0; i < 2; ++i)
      {
        for (int j = 0; j < 2; ++j)
        {
          read[i + 10*j + 0] = mass*velocity0[j];
          read[i + 10*j + 2] = mass*velocity1[j];
          read[i + 10*j + 4] = 0;
          read[i + 10*j + 6] = mass;
          read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                        + velocity1[j]*velocity1[j]);
          jacobian[j][2][2][i] = 1.;
        }

        jacobian[0][0][0][i] = 1.;
        jacobian[0][1][0][i] = 0.;
        jacobian[0][0][1][i] = 0.;
        jacobian[0][1][1][i] = 1.;

        jacobian[1][0][0][i] = -1.;
        jacobian[1][1][0][i] = 0.;
        jacobian[1][0][1][i] = 0.;
        jacobian[1][1][1][i] = -1.;
      }
      int i_axis_arg [] {1, 1};
      bool flip_arg [] {true, true};
      flip_arg[i_give] = false;
      cartdg::cpg_euler_hll_deformed<3, 2>(&read[0], &write[0],
                                           &(jacobian[0][0][0][0]), mult, i_axis_arg, flip_arg, 1.4);
      for (int j = 0; j < 2; ++j)
      {
        for (int i_var = 0; i_var < 5; ++i_var)
        {
          if (j == i_give)
          {
            REQUIRE(write[10*j + 2*i_var    ] == 0.);
            REQUIRE(write[10*j + 2*i_var + 1] == 0.);
          }
          else
          {
            double correct_d_flux = 3*340*read[2*i_var + 10*i_give] - 2*340*read[2*i_var + 10*(1 - i_give)];
            if (i_var == 4) correct_d_flux += 340*pressure;
            correct_d_flux *= 0.7;
            REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_d_flux).margin(1e-10));
            REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_d_flux).margin(1e-10));
          }
        }
      }
    }
  }

  SECTION("Oblique boundary easy")
  {
    for (int i_give : {0, 1})
    {
      double jacobian [2][3][3][2] {};
      double velocity0 [2] {};
      double velocity1 [2] {};
      velocity1[i_give] = 3*340*(1 - 2*i_give);
      velocity1[1 - i_give] = 2*340*(1 - 2*i_give);
      for (int i = 0; i < 2; ++i)
      {
        for (int j = 0; j < 2; ++j)
        {
          read[i + 10*j + 0] = mass*velocity0[j];
          read[i + 10*j + 2] = mass*velocity1[j];
          read[i + 10*j + 4] = 0;
          read[i + 10*j + 6] = mass;
          read[i + 10*j + 8] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                        + velocity1[j]*velocity1[j]);
          jacobian[j][2][2][i] = 1.;
        }

        jacobian[0][0][0][i] = 2.;
        jacobian[0][1][0][i] = 1.;
        jacobian[0][0][1][i] = 0;
        jacobian[0][1][1][i] = 1.;

        jacobian[1][0][0][i] = 2.;
        jacobian[1][1][0][i] = 1;
        jacobian[1][0][1][i] = 0.;
        jacobian[1][1][1][i] = 1.;
      }
      int i_axis_arg [] {1, 1};
      cartdg::cpg_euler_hll_deformed<3, 2>(&read[0], &write[0],
                                           &(jacobian[0][0][0][0]), mult, i_axis_arg, flip, 1.4);
      for (int j = 0; j < 2; ++j)
      {
        for (int i_var = 0; i_var < 5; ++i_var)
        {
          if (j == i_give)
          {
            REQUIRE(write[10*j + 2*i_var    ] == 0.);
            REQUIRE(write[10*j + 2*i_var + 1] == 0.);
          }
          else
          {
            double correct_d_flux = 3*340*read[2*i_var + 10*i_give] - 2*340*read[2*i_var + 10*(1 - i_give)];
            if (i_var == 4) correct_d_flux += 340*pressure;
            correct_d_flux *= 0.7;
            REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_d_flux).margin(1e-10));
            REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_d_flux).margin(1e-10));
          }
        }
      }
    }
  }

  SECTION("Oblique boundary difficult")
  {
    for (int i_give : {0, 1})
    {
      double jacobian [2][2][2][2] {};
      double velocity0 [2] {};
      double velocity1 [2] {};
      velocity0[i_give] = 1.5*340*(1 - 2*i_give);
      velocity1[i_give] = 1.5*680*(1 - 2*i_give);
      velocity0[1 - i_give] = 340*(1 - 2*i_give);
      velocity1[1 - i_give] = 680*(1 - 2*i_give);
      for (int i = 0; i < 2; ++i)
      {
        for (int j = 0; j < 2; ++j)
        {
          read[i + 8*j + 0] = mass*velocity0[j];
          read[i + 8*j + 2] = mass*velocity1[j];
          read[i + 8*j + 4] = mass;
          read[i + 8*j + 6] = pressure/0.4 + 0.5*mass*(  velocity0[j]*velocity0[j]
                                                       + velocity1[j]*velocity1[j]);
        }

        jacobian[0][0][0][i] = 2.;
        jacobian[0][1][0][i] = 1.;
        jacobian[0][0][1][i] = 0;
        jacobian[0][1][1][i] = 1.;

        jacobian[1][0][0][i] = 1.;
        jacobian[1][1][0][i] = 2.5;
        jacobian[1][0][1][i] = -2.;
        jacobian[1][1][1][i] = -1;
      }
      int i_axis_arg [] {1, 0};
      cartdg::cpg_euler_hll_deformed<2, 2>(&read[0], &write[0],
                                           &(jacobian[0][0][0][0]), mult, i_axis_arg, flip, 1.4);
      for (int j = 0; j < 2; ++j)
      {
        for (int i_var = 0; i_var < 4; ++i_var)
        {
          if (j == i_give)
          {
            REQUIRE(write[8*j + 2*i_var    ] == 0.);
            REQUIRE(write[8*j + 2*i_var + 1] == 0.);
          }
          else
          {
            double correct_d_flux = 1.5*3*340*read[2*i_var + 8*i_give] - 3*340*read[2*i_var + 8*(1 - i_give)];
            if (i_var == 3) correct_d_flux += 0.5*3*340*pressure;
            correct_d_flux *= 0.7;
            REQUIRE(write[8*j + 2*i_var    ] == Approx(correct_d_flux/(2*(1 + j))).margin(1e-10));
            REQUIRE(write[8*j + 2*i_var + 1] == Approx(correct_d_flux/(2*(1 + j))).margin(1e-10));
          }
        }
      }
    }
  }

  SECTION("Opposing supersonic flows")
  {
    double jacobian [2][3][3][2] {};
    double velocity0 [] {680, -680};
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        read[i + 10*j + 0] = mass*velocity0[j];
        read[i + 10*j + 2] = 0;
        read[i + 10*j + 4] = 0;
        read[i + 10*j + 6] = mass;
        read[i + 10*j + 8] = energy;

        for (int k = 0; k < 3; ++k)
        {
          jacobian[j][k][k][i] = 1.;
        }
      }
    }

    for (int i = 0; i < 20; ++i) write[i] = 0;
    int i_axis_arg [] {0, 0};
    cartdg::cpg_euler_hll_deformed<3, 2>(&read[0], &write[0],
                                         &(jacobian[0][0][0][0]), mult, i_axis_arg, flip, 1.4);
    REQUIRE(write[2*3     ] == Approx(0.7*mass*680));
    REQUIRE(write[2*3 + 10] == Approx(0.7*mass*680));
  }
}
