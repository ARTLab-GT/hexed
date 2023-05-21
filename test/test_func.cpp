#include <catch2/catch_all.hpp>

#include <hexed/config.hpp>
#include <hexed/Element_func.hpp>
#include <hexed/Qpoint_func.hpp>
#include <hexed/Domain_func.hpp>
#include <hexed/Spacetime_func.hpp>
#include <hexed/Surface_func.hpp>
#include <hexed/Gauss_lobatto.hpp>
#include <hexed/Gauss_legendre.hpp>
#include <hexed/Deformed_element.hpp>

class Arbitrary_func : public hexed::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return n_dim;}
  virtual std::string variable_name(int n_dim, int i_var) const {return "arbitrary" + std::to_string(i_var);}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    auto result = pos;
    double mult [] {-2, 0.3, 4, 0., 0.01};
    for (int i_dim = 0; i_dim < int(pos.size()); ++i_dim) {
      result[i_dim] = pos[i_dim]*mult[i_dim] - 2*time;
    }
    return result;
  }
};

std::vector<std::vector<double>> test_pos
{
  {},
  {0.2, 7.9},
  {0.},
  {0., 0., 0.},
};

std::vector<double> test_time {1, 0., 10, -1.3};

std::vector<std::vector<double>> test_state
{
  {},
  {1.3, 0.02},
  {1.},
  {0., 0., 0.},
};

std::vector<std::vector<double>> test_normal
{
  {},
  {3., 4.},
  {2.},
  {0., 0., 1.},
};

std::vector<std::vector<double>> test_error
{
  {},
  {1.3 - (-2*0.2 - 2*0), 0.02 - (0.3*7.9 - 2*0)},
  {1. - (-2*0. - 2*10.)},
  {2*-1.3, 2*-1.3, 2*-1.3},
};

TEST_CASE("Constant_func")
{
  std::vector<double> value {0.3, -0.7};
  hexed::Constant_func cf (value);
  REQUIRE(cf.n_var(2) == 2);
  REQUIRE(cf.n_var(3) == 2);
  for (auto pos : test_pos) {
    for (auto time : test_time) {
      auto result = cf(pos, time);
      REQUIRE(result == value);
    }
  }
}

TEST_CASE("Error_func")
{
  Arbitrary_func af;
  hexed::Error_func ef(af);
  REQUIRE(ef.n_var(2) == 2);
  REQUIRE(ef.n_var(3) == 3);
  REQUIRE(ef.variable_name(2, 1) == "state1_errsq");
  for (unsigned i_test = 0; i_test < test_pos.size(); ++i_test) {
    auto error = ef(test_pos[i_test], test_time[i_test], test_state[i_test]);
    REQUIRE(error.size() == test_error[i_test].size());
    for (unsigned i_val = 0; i_val < error.size(); ++i_val) {
      REQUIRE(std::sqrt(error[i_val]) == std::abs(test_error[i_test][i_val]));
    }
  }
}

TEST_CASE("Inherited methods")
{
  SECTION("Qpoint_func")
  {
    Arbitrary_func af;
    hexed::Gauss_lobatto basis {2};
    hexed::Storage_params params {2, 4, 2, 2};
    hexed::Deformed_element elem {params, {0, 1}};
    hexed::Qpoint_func& af_ref {af};
    auto result = af_ref(elem, basis, 1, 5.);
    REQUIRE(result.size() == 2);
    REQUIRE(result[1] == Catch::Approx(0.3*2. - 10.));
  }
  SECTION("Surface_func")
  {
    // Verify that (Surface_func&)(State_variables&) gives you state variables
    // regardless of the other inputs
    hexed::State_variables sv;
    hexed::Surface_func& surface {sv};
    for (auto pos : test_pos) {
      for (auto time : test_time) {
        for (auto state : test_state) {
          for (auto normal : test_normal) {
            REQUIRE(surface(pos, time, state, normal) == state);
          }
        }
      }
    }
  }
}

TEST_CASE("Vortex")
{
  SECTION("Solves inviscid flow equations")
  {
    std::vector<double> freestream {9., -0.2, 0.8, 2e5};
    hexed::Isentropic_vortex vortex(freestream);
    REQUIRE(vortex.n_var(2) == 4);

    std::vector<double> test_pos []
    {
      {0.01, 0.05},
      {-0.03, -0.03},
      {0., 0.24},
      {0, -0.01},
      {0., 0.},
      {2, 10},
    };
    double test_time [] {0., 0., 0.11, -0.03, 0.2, 3};

    for (int i_test = 0; i_test < 5; ++i_test)
    {
      std::vector<double> flux_grad [2];
      flux_grad[0].resize(4); flux_grad[1].resize(4);
      std::vector<double> rate (4);
      double diff = 1.e-6;
      for (int dir : {-1, 1})
      {
        for (int i_dim : {0, 1})
        {
          auto pos = test_pos[i_test];
          pos[i_dim] += dir*diff;
          auto state = vortex(pos, test_time[i_test]);
          double pressure = state[3] - 0.5/state[2]*(state[0]*state[0] + state[1]*state[1]);
          pressure *= (vortex.heat_rat - 1.);
          std::vector<double> flux (4);
          for (int j_dim = 0; j_dim < 2; ++j_dim)
          {
            flux[j_dim] += state[i_dim]*state[j_dim]/state[2];
          }
          flux[i_dim] += pressure;
          flux[2] = state[i_dim];
          flux[3] = state[i_dim]/state[2]*(state[3] + pressure);
          for (int i_var = 0; i_var < 4; ++i_var)
          {
            flux_grad[i_dim][i_var] += dir*flux[i_var]/(2*diff);
          }
        }
        auto pos = test_pos[i_test];
        auto state = vortex(pos, test_time[i_test] + dir*diff);
        for (int i_var = 0; i_var < 4; ++i_var)
        {
          rate[i_var] += dir*state[i_var]/(2*diff);
        }
      }

      for (int i_var = 0; i_var < 4; ++i_var)
      {
        REQUIRE(rate[i_var] + flux_grad[0][i_var] + flux_grad[1][i_var]
                == Catch::Approx(0.).margin(1.e-3*std::abs(freestream[i_var])));
      }
    }
  }

  SECTION("Max tangential velocity is correct")
  {
    std::vector<double> freestream = {0., 0., 1., 2e5};
    hexed::Isentropic_vortex vortex (freestream);
    auto state = vortex(std::vector<double> {0., vortex.argmax_radius}, 0.);
    double sound_speed = std::sqrt(vortex.heat_rat*(vortex.heat_rat - 1.)*freestream[3]);
    double tang_veloc = -state[0]/state[2];
    REQUIRE(tang_veloc/sound_speed == Catch::Approx(vortex.max_nondim_veloc));
  }
}

TEST_CASE("Doublet")
{
  double veloc [] {10., -1.};
  double angle_of_attack {std::atan2(veloc[1], veloc[0])};
  double speed_sq {veloc[0]*veloc[0] + veloc[1]*veloc[1]};
  double mass {1.2};
  double pres {101000.};
  double ener {pres/0.3 + 0.5*mass*speed_sq};
  std::vector<double> freestream = {veloc[0]*mass, veloc[1]*mass, mass, ener};
  hexed::Doublet doublet {freestream};
  REQUIRE(doublet.n_var(2) == 4);
  doublet.radius = 0.8;
  doublet.heat_rat = 1.3; // just to make sure 1.4 isn't hardcoded anywhere
  doublet.location = {0.05, 0.03};
  for (double angle_degrees : {0., 10., 90., 160.})
  {
    // compute state on circle for some arbitrary time (time doesn't matter).
    double angle {angle_degrees*M_PI/180.};
    double sin {std::sin(angle)};
    double cos {std::cos(angle)};
    std::vector<double> state {doublet({doublet.radius*cos + 0.05, doublet.radius*sin + 0.03}, 0.2089)};
    double momentum_magnitude = std::sqrt(state[0]*state[0] + state[1]*state[1]);
    // require velocity is tangential to surface
    REQUIRE(state[0]*cos + state[1]*sin == Catch::Approx(0.).scale(10.));
    // require pressure approximately satisfies Bernoulli eq.
    double pres_coef = (0.3*(state[3] - 0.5*momentum_magnitude*momentum_magnitude/state[2])
                        - pres)/(0.5*mass*speed_sq);
    REQUIRE(pres_coef == Catch::Approx(1 - 4*std::pow(std::sin(angle - angle_of_attack), 2)).epsilon(1e-2));
    // require state at large radius is approximately freestream
    state = doublet({1e3*doublet.radius*cos + 0.05, 1e3*doublet.radius*sin + 0.03}, 0.2089);
    REQUIRE(state[0] == Catch::Approx(mass*veloc[0]).epsilon(1e-3));
    REQUIRE(state[1] == Catch::Approx(mass*veloc[1]).epsilon(1e-3));
    REQUIRE(state[2] == Catch::Approx(mass         ).epsilon(1e-3));
    REQUIRE(state[3] == Catch::Approx(ener         ).epsilon(1e-3));
  }
}

TEST_CASE("Pressure_stress")
{
  std::vector<double> pressure {1e5, 2e4, 3.2e4};
  std::vector<double> pos {}; // size of `pos` shouldn't matter
  double time {0.2};
  double veloc [3] {0.2, -2., 10.};
  double mass {0.7};
  std::vector<std::vector<double>> unit_normal {
    {},
    {3./5., 4./5.},
    {1.},
    {0., 0., 1.},
  };
  hexed::Pressure_stress fpa {1.2};
  // verify that when you back out the state,
  // `Pressure_stress` gives you pressure times unit normal
  for (unsigned i_normal = 0; i_normal < test_normal.size(); ++i_normal)
  {
    auto normal {test_normal[i_normal]};
    for (double pres : pressure)
    {
      std::vector<double> state;
      double kin_ener = 0.;
      for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
        kin_ener += 0.5*mass*veloc[i_dim]*veloc[i_dim];
        state.push_back(mass*veloc[i_dim]);
      }
      state.push_back(mass);
      state.push_back(pres/0.2 + kin_ener);
      auto computed {fpa(pos, time, state, normal)};
      REQUIRE(computed.size() == normal.size());
      for (unsigned i_dim = 0; i_dim < normal.size(); ++i_dim) {
        REQUIRE(-computed[i_dim]/pres == Catch::Approx(unit_normal[i_normal][i_dim]).scale(1.));
      }
    }
  }
}

class Nonsmooth_test : public hexed::Spacetime_func
{
  public:
  int n_var(int n_dim) const {return 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n = 10000;
    return {double(std::rand()%n)/n, std::exp(pos[1]/.3)};
  }
};

TEST_CASE("elementwise functions")
{
  hexed::Gauss_legendre basis(hexed::config::max_row_size);
  hexed::Deformed_element elem({2, 4, 2, hexed::config::max_row_size}, {}, .3);
  elem.vertex(1).pos[1] *= 2.;
  elem.vertex(3).pos[1] *= 2.;
  double faces [6][6*hexed::config::max_row_size];
  for (int i_face = 0; i_face < 6; ++i_face) elem.faces[i_face] = faces[i_face];
  elem.set_jacobian(basis);
  Arbitrary_func func;
  SECTION("average")
  {
    hexed::Elem_average avg(func);
    auto result = avg(elem, basis, 7.);
    REQUIRE(result.size() == 2);
    REQUIRE(result[0] == Catch::Approx(-2*.3/2. - 2*7.));
    REQUIRE(result[1] == Catch::Approx(.3*.6/2. - 2*7.));
  }
  SECTION("L2 norm")
  {
    hexed::Elem_l2 l2(func);
    auto result = l2(elem, basis, 7.);
    REQUIRE(result.size() == 2);
    REQUIRE(result[0] == Catch::Approx(std::sqrt(2.*2.*.3*.3/3. + 2*2.*.3/2.*2*7. + 14*14)));
    REQUIRE(result[1] == Catch::Approx(std::sqrt(.3*.3*.6*.6/3. - 2*.3*.6/2.*2*7.+ 14*14)));
  }
  SECTION("nonsmoothness")
  {
    srand(406);
    Nonsmooth_test nonsm;
    auto result = hexed::Elem_nonsmooth(nonsm)(elem, basis, .7);
    auto norm = hexed::Elem_l2(nonsm)(elem, basis, .7);
    REQUIRE(result.size() == 2);
    REQUIRE(norm.size() == 2);
    // results should not be greater than the l2 norm
    for (int i = 0; i < 2; ++i) {
      result[i] /= norm[i];
      REQUIRE(result[i] <= 1 + 1e-14);
    }
    REQUIRE(result[0] > .1); // random input should have nonsmoothness on the order of 1.
    REQUIRE(result[1] < .01); // analytic input should have small nonsmoothness
  }
}

TEST_CASE("Spatial_gaussian")
{
  REQUIRE_THROWS(hexed::Spatial_gaussian({}));
  hexed::Spatial_gaussian gaus({.8, .5});
  REQUIRE(gaus({0., 0.}, 0.).size() == 1);
  REQUIRE(gaus({0., 0.}, 1.7)[0] == Catch::Approx(1.));
  REQUIRE(gaus({.8, 0., 0.}, 1.7)[0]/gaus({0., .5, 0.}, 1.)[0] == Catch::Approx(1.));
  REQUIRE(gaus({.8, 0., 0.}, 1.7)[0]/gaus({0., 0., .5}, 1.)[0] == Catch::Approx(1.));
  REQUIRE(gaus({0., .5, 0.}, 1.7)[0]/gaus({0., 1., 0.}, 1.)[0] == Catch::Approx(std::exp(1.5)));
}

TEST_CASE("Pow")
{
  hexed::Gauss_lobatto basis(2);
  hexed::Element elem({2, 4, 2, 2}, {}, 10. + 2./3.);
  Arbitrary_func func;
  hexed::Pow p(func, 3);
  REQUIRE(p.n_var(2) == 2);
  REQUIRE(p.variable_name(2, 1) == "(arbitrary1)^3");
  auto result = p(elem, basis, 1, .1);
  REQUIRE(result.size() == 2);
  REQUIRE(result[0] == Catch::Approx(-.008));
  REQUIRE(result[1] == Catch::Approx(27));
}

TEST_CASE("Qf_concat")
{
  Arbitrary_func af;
  hexed::Constant_func c({1., 2.});
  hexed::Qf_concat concat({&af, &c});
  // check that number of variables is sum of `af` and `c`
  REQUIRE(concat.n_var(1) == 3);
  REQUIRE(concat.n_var(3) == 5);
  // check that names of variables of `af` and `c` are retained
  REQUIRE(concat.variable_name(3, 0) == "arbitrary0");
  REQUIRE(concat.variable_name(3, 2) == "arbitrary2");
  REQUIRE(concat.variable_name(3, 3) == "constant0");
  // check that values are correct
  hexed::Gauss_lobatto basis {2};
  hexed::Storage_params params {2, 4, 2, 2};
  hexed::Deformed_element elem {params, {0, 1}};
  auto result = concat(elem, basis, 1, 5.);
  REQUIRE(result.size() == 4);
  REQUIRE(result[1] == Catch::Approx(0.3*2. - 10.));
  REQUIRE(result[2] == Catch::Approx(1.));
  REQUIRE(result[3] == Catch::Approx(2.));
}

TEST_CASE("Mach")
{
  std::vector<double> state {3.*1.2, 4.*1.2, 1.2, 101325/.4 + .5*25*1.2};
  auto m = hexed::Mach()(std::vector<double>{1., 2.}, 0., state);
  REQUIRE(m.size() == 1);
  REQUIRE(m[0] == Catch::Approx(5./std::sqrt(1.4*101325/1.2)));
}

TEST_CASE("Equiangle_skewness")
{
  hexed::Gauss_legendre basis(2);
  hexed::Deformed_element elem({2, 4, 2, 2});
  double faces[4][2*2*4];
  for (int i_face = 0; i_face < 4; ++i_face) elem.faces[i_face] = faces[i_face];
  hexed::Equiangle_skewness skew;
  elem.set_jacobian(basis);
  REQUIRE(skew(elem, basis, 0.)[0] == Catch::Approx(0.).scale(1.));
  elem.vertex(1).pos = {1./std::sqrt(3.), 1., 0.};
  elem.set_jacobian(basis);
  REQUIRE(skew(elem, basis, 0.)[0] == Catch::Approx(1./3.));
  elem.vertex(1).pos = {.5, .5, 0.};
  elem.set_jacobian(basis);
  REQUIRE(skew(elem, basis, 0.)[0] == Catch::Approx(1.));
}

#if HEXED_USE_NLOPT
TEST_CASE("Ringleb")
{
  hexed::Ringleb ringleb;
  std::vector<double> pos {1.13206394e+01, 4.02855640e+00};
  auto state = ringleb(pos, 0.);
  REQUIRE(std::sqrt(state[0]*state[0] + state[1]*state[1])/state[2] == Catch::Approx(.2));
};
#endif
