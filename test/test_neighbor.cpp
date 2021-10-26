#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <neighbor/read_copy.hpp>
#include <neighbor/write_copy.hpp>
#include <neighbor/hll_cpg_euler.hpp>
#include <neighbor/ausm_plus_up_cpg_euler.hpp>
#include <neighbor/variable_jump.hpp>
#include <get_neighbor_derivative.hpp>
#include <get_cont_visc.hpp>
#include <Storage_params.hpp>
#include <Gauss_lobatto.hpp>

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

TEST_CASE("neighbor kernel write_copy<>()")
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

TEST_CASE("hll_cpg_euler")
{
  double mass = 1.225;
  double pressure = 101325;
  double read[20] {};
  double write[20] {};
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
    cartdg::hll_cpg_euler<3, 2>(&read[0], &write[0], 1., 0, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_flux = read[2*i_var];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= 3*340;
        if (i_var == 0) correct_flux += pressure;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) write[i] = 0;
    cartdg::hll_cpg_euler<3, 2>(&read[0], &write[0], 1., 1, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_flux = read[2*i_var + 10];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= -3*340;
        if (i_var == 1) correct_flux += pressure;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
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
    cartdg::hll_cpg_euler<3, 2>(&read[0], &write[0], 1., 0, 1.4);
    REQUIRE(write[2*3     ] == Approx(0).margin(1e-8));
    REQUIRE(write[2*3 + 10] == Approx(0).margin(1e-8));
  }
}

TEST_CASE("ausm_plus_up_cpg_euler")
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
    cartdg::ausm_plus_up_cpg_euler<3, 2>(&read[0], &write[0], mult, 0, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_flux = read[2*i_var];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= 3*340;
        if (i_var == 0) correct_flux += pressure;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
      }
    }
    for (int i = 0; i < 20; ++i) write[i] = 0;
    cartdg::ausm_plus_up_cpg_euler<3, 2>(&read[0], &write[0], mult, 1, 1.4);
    for (int j = 0; j < 2; ++j)
    {
      for (int i_var = 0; i_var < 5; ++i_var)
      {
        double correct_flux = read[2*i_var + 10];
        if (i_var == 4) correct_flux += pressure;
        correct_flux *= -3*340;
        if (i_var == 1) correct_flux += pressure;
        REQUIRE(write[10*j + 2*i_var    ] == Approx(correct_flux).margin(1.e-8));
        REQUIRE(write[10*j + 2*i_var + 1] == Approx(correct_flux).margin(1.e-8));
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
    cartdg::ausm_plus_up_cpg_euler<3, 2>(&read[0], &write[0], mult, 0, 1.4);
    REQUIRE(write[2*3     ] == Approx(0).margin(1e-8));
    REQUIRE(write[2*3 + 10] == Approx(0).margin(1e-8));
  }
}

TEST_CASE("jump kernel")
{
  const int row_size = CARTDG_MAX_BASIS_ROW_SIZE;
  const int n_qpoint = row_size*row_size*row_size;
  double read  [3][1][row_size][row_size][row_size];
  double write [3][1][row_size][row_size][row_size];
  double weight = 0.2;
  for (int i = 0; i < n_qpoint; ++i)
  {
    *(&read[0][0][0][0][0] + i) = 0.5;
    *(&read[1][0][0][0][0] + i) = 2.;
    *(&read[2][0][0][0][0] + i) = 0.7;
    for (int j = 0; j < 3; ++j)
    {
      *(&write[j][0][0][0][0] + i) = 0.1;
    }
  }

  cartdg::variable_jump<n_qpoint, row_size>({&read[2][0][0][0][0], &read[1][0][0][0][0]}, {&write[1][0][0][0][0], &write[0][0][0][0][0]}, 1, weight);
  cartdg::variable_jump<n_qpoint, row_size>({&read[0][0][0][0][0], &read[2][0][0][0][0]}, {&write[2][0][0][0][0], &write[1][0][0][0][0]}, 1, weight);
  double correct = (0.7 - 0.5)/2./0.2 + 0.1;
  REQUIRE(write[2][0][0][row_size - 1][0] == Approx(correct));
  REQUIRE(write[2][0][row_size - 1][row_size - 1][row_size - 1] == Approx(correct));
  REQUIRE(write[1][0][0][0][0] == Approx(correct));
  REQUIRE(write[1][0][row_size - 1][0][row_size - 1] == Approx(correct));
  REQUIRE(write[0][0][0][row_size - 1][0] == 0.1);
  REQUIRE(write[2][0][0][0][0] == 0.1);
  REQUIRE(write[1][0][1][1][1] == 0.1); // might fail if CARTDG_MAX_BASIS_ROW_SIZE == 2
}

TEST_CASE("neighbor_derivative")
{
  const int row_size = CARTDG_MAX_BASIS_ROW_SIZE;
  const int n_qpoint = row_size*row_size*row_size;
  cartdg::Storage_params params {1, 5, 3, row_size};
  cartdg::Kernel_settings settings;
  cartdg::Gauss_lobatto basis {row_size};
  cartdg::elem_vec elements;
  for (int i_elem = 0; i_elem < 3; ++i_elem)
  {
    elements.emplace_back(new cartdg::Element {params});
  }
  for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
  {
    elements[0]->stage(0)[2*n_qpoint + i_qpoint] = 0.5;
    elements[1]->stage(0)[2*n_qpoint + i_qpoint] = 2.;
    elements[2]->stage(0)[2*n_qpoint + i_qpoint] = 0.7;
    for (int i_elem = 0; i_elem < 3; ++i_elem)
    {
      elements[i_elem]->derivative()[i_qpoint] = 0.1;
    }
  }
  cartdg::elem_con_vec connections {{}, {{elements[2].get(), elements[1].get()}, {elements[0].get(), elements[2].get()}}, {}};
  cartdg::get_neighbor_derivative(3, row_size)(connections, 2, 1, basis, settings);
  double correct = (0.7 - 0.5)/2./basis.node_weights()[0] + 0.1;
  int end = (row_size - 1)*(row_size*row_size + 1);
  REQUIRE(elements[0]->derivative()[(row_size - 1)*row_size + 0  ] == Approx(correct));
  REQUIRE(elements[0]->derivative()[(row_size - 1)*row_size + end] == Approx(correct));
  REQUIRE(elements[2]->derivative()[(0           )*row_size + 0  ] == Approx(correct));
  REQUIRE(elements[2]->derivative()[(0           )*row_size + end] == Approx(correct));
  REQUIRE(elements[1]->derivative()[(row_size - 1)*row_size + 0 ] == 0.1);
  REQUIRE(elements[0]->derivative()[(0           )*row_size + 0 ] == 0.1);
  REQUIRE(elements[2]->derivative()[(1           )*row_size + row_size*row_size + 1] == 0.1); // might fail if CARTDG_MAX_BASIS_ROW_SIZE == 2
}

TEST_CASE("continuous viscosity kernel")
{
  cartdg::Kernel_settings settings;
  cartdg::Storage_params params {1, 1, 2, 4};
  cartdg::elem_vec elements;
  for (int i = 0; i < 4; ++i) elements.emplace_back(new cartdg::Element {params});
  /*
  element order is transposed, just for fun:
    21 23   31 33
    20 22   30 32
  ^ 
  | 01 03   11 13
  y 00 02   10 12
   x --->
  */
  cartdg::elem_con_vec connections {{{elements[0].get(), elements[1].get()},
                                     {elements[2].get(), elements[3].get()}},
                                    {{elements[0].get(), elements[2].get()},
                                     {elements[1].get(), elements[3].get()}}};
  for (int i_point = 0; i_point < 4; ++i_point)
  {
    elements[0]->viscosity()[i_point] = 0.;
    elements[1]->viscosity()[i_point] = 0.5;
    elements[2]->viscosity()[i_point] = 2.;
    elements[3]->viscosity()[i_point] = 1.;
  }

  cartdg::get_cont_visc(2, CARTDG_MAX_BASIS_ROW_SIZE)(connections, settings);
  REQUIRE(elements[0]->viscosity()[0] == Approx(0.0));
  REQUIRE(elements[0]->viscosity()[1] == Approx(2.0));
  REQUIRE(elements[0]->viscosity()[2] == Approx(0.5));
  REQUIRE(elements[0]->viscosity()[3] == Approx(2.0));
  REQUIRE(elements[1]->viscosity()[0] == Approx(0.5));
  REQUIRE(elements[1]->viscosity()[1] == Approx(2.0));
  REQUIRE(elements[1]->viscosity()[2] == Approx(0.5));
  REQUIRE(elements[1]->viscosity()[3] == Approx(1.0));
  REQUIRE(elements[2]->viscosity()[0] == Approx(2.0));
  REQUIRE(elements[2]->viscosity()[1] == Approx(2.0));
  REQUIRE(elements[3]->viscosity()[0] == Approx(2.0));
  REQUIRE(elements[3]->viscosity()[1] == Approx(2.0));
  REQUIRE(elements[3]->viscosity()[2] == Approx(1.0));
  REQUIRE(elements[3]->viscosity()[3] == Approx(1.0));
}
