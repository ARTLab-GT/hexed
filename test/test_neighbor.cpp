#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <neighbor/hll_cpg_euler.hpp>
#include <neighbor/ausm_plus_up_cpg_euler.hpp>
#include <get_nonpen_convective.hpp>
#include <get_prolong.hpp>
#include <get_restrict.hpp>
#include <Storage_params.hpp>
#include <Gauss_lobatto.hpp>
#include <Gauss_legendre.hpp>
#include <Kernel_settings.hpp>

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

TEST_CASE("nonpen")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Gauss_legendre basis {row_size};
  cartdg::Kernel_settings settings;
  cartdg::Storage_params params {1, 4, 2, row_size};
  cartdg::Deformed_element elem {params, {0, 0, 0}};
  cartdg::Face_index index {&elem, 0, 1};
  cartdg::Deformed_elem_wall wall {index, 0};
  for (int i_qpoint = 0; i_qpoint < row_size; ++i_qpoint) {
    double qpoint_jacobian [] {1., 3., 0., 4.};
    for (int i_jac = 0; i_jac < 4; ++i_jac) wall.jacobian()[i_jac*row_size + i_qpoint] = qpoint_jacobian[i_jac];
    double state [] {1., 1., 1.2, 1e5/0.4 + 0.5*1.2*2.};
    for (int i_var = 0; i_var < 4; ++i_var) elem.face()[(4 + i_var)*row_size + i_qpoint] = state[i_var];
  }
  elem.face()[4*row_size] = 0.; // momentum[0] at qpoint[0]
  elem.face()[5*row_size] = 0.; // momentum[1] at qpoint[0]
  cartdg::def_elem_wall_vec vec {wall};
  cartdg::get_nonpen_convective(2, row_size)(vec, basis, settings);
  double normal_momentum_flux = 4./5.*elem.face()[4*row_size] - 3./5.*elem.face()[5*row_size];
  double tangential_momentum_flux = 3./5.*elem.face()[4*row_size] + 4./5.*elem.face()[5*row_size];
  REQUIRE(normal_momentum_flux == Approx(1e5*5.)); // when momentum == 0, normal momentum flux is pressure
  REQUIRE(tangential_momentum_flux == Approx(0.).scale(1.)); // when momentum != 0, tangential momentum flux is 0
  REQUIRE(elem.face()[6*row_size + 1] == Approx(0.).scale(1.)); // when momentum != 0, mass flux is 0
}

TEST_CASE("prolong/restrict")
{
  // test that prolongation/restriction operators are approximately correct for an exponential function
  const int row_size {cartdg::config::max_row_size};
  cartdg::Kernel_settings settings;
  cartdg::Storage_params params {2, 5, 3, row_size};
  cartdg::Gauss_legendre basis {row_size};

  double coarse [5][row_size][row_size] {};
  cartdg::ref_face_vec ref_faces;
  for (int i_dim = 0; i_dim < 3; ++i_dim) ref_faces.push_back({});
  ref_faces[2].emplace_back(new cartdg::Refined_face {params, coarse[0][0]});

  SECTION("prolong")
  {
    for (int i_var = 0; i_var < 5; ++i_var) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          coarse[i_var][i_node][j_node] = std::exp(basis.node(i_node) + 0.5*basis.node(j_node)) + i_var;
        }
      }
    }
    cartdg::get_prolong(3, row_size)(ref_faces, basis, settings);
    for (int i_half : {0, 1}) {
      for (int j_half : {0, 1}) {
        for (int i_node = 0; i_node < row_size; ++i_node) {
          for (int j_node = 0; j_node < row_size; ++j_node) {
            for (int i_var = 0; i_var < 5; ++i_var) {
              double prolonged {ref_faces[2][0]->fine_face(i_half*2 + j_half)[(i_var*row_size + i_node)*row_size + j_node]};
              double correct {std::exp((basis.node(i_node) + i_half)/2. + 0.5*(basis.node(j_node) + j_half)/2.) + i_var};
              REQUIRE(prolonged == Approx(correct).margin(1e-4));
            }
          }
        }
      }
    }
  }

  SECTION("restrict")
  {
    for (int i_half : {0, 1}) {
      for (int j_half : {0, 1}) {
        for (int i_node = 0; i_node < row_size; ++i_node) {
          for (int j_node = 0; j_node < row_size; ++j_node) {
            for (int i_var = 0; i_var < 5; ++i_var) {
              ref_faces[2][0]->fine_face(i_half*2 + j_half)[(i_var*row_size + i_node)*row_size + j_node]
                = std::exp((basis.node(i_node) + i_half)/2. + 0.5*(basis.node(j_node) + j_half)/2.) + i_var;
            }
          }
        }
      }
    }
    cartdg::get_restrict(3, row_size)(ref_faces, basis, settings);
    for (int i_var = 0; i_var < 5; ++i_var) {
      for (int i_node = 0; i_node < row_size; ++i_node) {
        for (int j_node = 0; j_node < row_size; ++j_node) {
          double restricted {coarse[i_var][i_node][j_node]};
          double correct {std::exp(basis.node(i_node) + 0.5*basis.node(j_node)) + i_var};
          REQUIRE(restricted == Approx(correct).margin(1e-4));
        }
      }
    }
  }
}
