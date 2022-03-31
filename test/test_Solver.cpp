#include <catch2/catch.hpp>
#include <cartdgConfig.hpp>
#include <Solver.hpp>

class Arbitrary_initializer : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 4;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {pos[0]*pos[1], 1., 2., 3.};
  }
};

class Bad_initializer : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {0., 0.};
  }
};

class Arbitrary_integrand : public cartdg::Domain_func
{
  public:
  virtual int n_var(int n_dim) const {return 3;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state) const
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0.};
  }
};

constexpr double velocs [3] {0.3, -0.7, 0.8};
constexpr double wave_number [3] {-0.1, 0.3, 0.2};
double scaled_veloc(int n_dim)
{
  double result = 0.;
  for (int i_dim = 0; i_dim < n_dim; ++i_dim) result += velocs[i_dim]*wave_number[i_dim];
  return result;
}

class Nonuniform_mass : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return n_dim + 2;};

  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n_dim {int(pos.size())};
    std::vector<double> result;
    double scaled_pos = 0.;
    double veloc_mag_sq = 0.;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      scaled_pos += wave_number[i_dim]*pos[i_dim];
      veloc_mag_sq += velocs[i_dim]*velocs[i_dim];
    }
    double mass = 1. + 0.1*std::sin(scaled_pos);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) result.push_back(mass*velocs[i_dim]);
    result.push_back(mass);
    result.push_back(1e5/0.4 + mass*0.5*veloc_mag_sq);
    return result;
  }
};

class Nonuniform_residual : public cartdg::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return n_dim + 2;};

  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n_dim {int(pos.size())};
    std::vector<double> result;
    double scaled_pos = 0.;
    double veloc_mag_sq = 0.;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) {
      scaled_pos += wave_number[i_dim]*pos[i_dim];
      veloc_mag_sq += velocs[i_dim]*velocs[i_dim];
    }
    double d_mass = -0.1*scaled_veloc(n_dim)*std::cos(scaled_pos);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) result.push_back(d_mass*velocs[i_dim]);
    result.push_back(d_mass);
    result.push_back(d_mass*0.5*veloc_mag_sq);
    return result;
  }
};

TEST_CASE("Solver")
{
  static_assert (cartdg::config::max_row_size >= 3); // this test was written for row size 3
  cartdg::Solver sol {2, 3, 0.8};

  SECTION("initialization and sampling")
  {
    sol.mesh().add_element(0, false, {0, 0, 0});
    int sn0 = sol.mesh().add_element(0, false, {1, 0, 0});
    int sn1 = sol.mesh().add_element(2, true, {-1, 0, 0});
    REQUIRE_THROWS(sol.initialize(Bad_initializer())); // if number of variables of func is wrong, should throw
    sol.initialize(Arbitrary_initializer());
    auto sample = sol.sample(0, false, sn0, 4, cartdg::State_variables()); // sample the midpoint of the element because we know the exact position
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Approx(1.2*0.4));
    REQUIRE(sample[1] == Approx(1.));
    REQUIRE(sample[2] == Approx(2.));
    REQUIRE(sample[3] == Approx(3.));
    sample = sol.sample(2, true, sn1, 4, cartdg::State_variables());
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Approx(-0.1*0.1));
  }

  SECTION("vertex relaxation")
  {
    int sn0 = sol.mesh().add_element(0, false, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0,  true, {1, 0, 0});
    sol.mesh().connect_cartesian(0, {sn0, sn1}, {0}, {false, true});
    sol.relax_vertices();
    REQUIRE(sol.sample(0, false, sn0, 4, cartdg::Position_func())[0] == Approx(0.8*0.5));
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Position_func())[0] == Approx(0.8*1.375));
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Position_func())[1] == Approx(0.8*0.5));
  }

  SECTION("local time step scale")
  {
    int sn0 = sol.mesh().add_element(0,  true, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0,  true, {0, 0, 0});
    int sn2 = sol.mesh().add_element(0, false, {-1, 0, 0});
    int sn3 = sol.mesh().add_element(1, false, {-2, -1, 0});
    int sn4 = sol.mesh().add_element(1, false, {-1, -1, 0});
    sol.mesh().connect_deformed(0, {sn0, sn1}, {{0, 0}, {1, 0}});
    sol.mesh().connect_cartesian(0, {sn2, sn0}, {0}, {false, true});
    sol.mesh().connect_hanging_cartesian(0, sn2, {sn3, sn4}, {1}, 0);
    sol.calc_jacobian();
    sol.set_local_tss();
    // in sn1, TSS is 0.5 because the element is stretched by a factor of 0.5
    REQUIRE(sol.sample(0,  true, sn1, 4, cartdg::Time_step_scale_func())[0] == Approx(0.5));
    // in sn2, TSS varies linearly between 1 and 0.5 (because it must be continuous with element sn1)
    REQUIRE(sol.sample(0, false, sn2, 4, cartdg::Time_step_scale_func())[0] == Approx(0.75));
    // TSS at hanging nodes should be set to match coarse element
    REQUIRE(sol.sample(1, false, sn3, 4, cartdg::Time_step_scale_func())[0] == Approx(1. - 0.0625));
    REQUIRE(sol.sample(1, false, sn4, 4, cartdg::Time_step_scale_func())[0] == Approx(0.5*(1. + 0.625)));
  }

  SECTION("field integrals")
  {
    SECTION("simple function, complex mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.mesh().add_element(0, false, {1, 0, 0});
      sol.mesh().add_element(1, false, {2, 2, 0});
      int sn0 = sol.mesh().add_element(0, true, {2, 0, 0});
      int sn1 = sol.mesh().add_element(0, true, {4, 0, 0});
      // connecting deformed elements with an empty space in between stretches them by a factor of 1.5
      sol.mesh().connect_deformed(0, {sn0, sn1}, {{0, 0}, {1, 0}});
      sol.calc_jacobian();
      std::vector<double> state {0.3, -10., 0.7, 32.};
      sol.initialize(cartdg::Constant_func(state));
      auto integral = sol.integral_field(cartdg::State_variables());
      REQUIRE(integral.size() == 4);
      double area = 5.25*0.8*0.8;
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Approx(state[i_var]*area));
      }
    }
    SECTION("complex function, simple mesh")
    {
      sol.mesh().add_element(0, false, {0, 0, 0});
      sol.initialize(cartdg::Constant_func({0.3, 0., 0., 0.}));
      auto integral = sol.integral_field(Arbitrary_integrand());
      REQUIRE(integral.size() == 3);
      REQUIRE(integral[0] == Approx(std::pow(0.8, 3)/3*std::pow(0.8, 4)/4 - 0.8*0.8*0.3));
    }
  }
}

class Test_mesh
{
  public:
  struct elem_handle {int ref_level; bool is_deformed; int serial_n;};
  virtual cartdg::Solver& solver() = 0;
  virtual std::vector<elem_handle> construct(cartdg::Boundary_condition* bc) = 0;
};

// creates a 2x2x2 mesh
class All_cartesian : public Test_mesh
{
  cartdg::Solver sol;
  public:

  All_cartesian()
  : sol{3, cartdg::config::max_row_size, 1.}
  {}

  virtual cartdg::Solver& solver() {return sol;}

  virtual std::vector<elem_handle> construct(cartdg::Boundary_condition* bc)
  {
    std::vector<elem_handle> handles;
    int bc_sn = sol.mesh().add_boundary_condition(bc);
    for (int i_elem = 0; i_elem < 8; ++i_elem) {
      std::vector<int> strides {4, 2, 1};
      std::vector<int> inds;
      for (int i_dim = 0; i_dim < 3; ++i_dim) inds.push_back((i_elem/strides[i_dim])%2);
      int sn = sol.mesh().add_element(0, false, inds);
      handles.push_back({0, false, sn});
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        if (inds[i_dim]) sol.mesh().connect_cartesian(0, {handles[i_elem - strides[i_dim]].serial_n, sn}, {i_dim});
        sol.mesh().connect_boundary(0, false, sn, i_dim, inds[i_dim], bc_sn);
      }
    }
    return handles;
  }
};

// creates a 3x3 mesh with vertex perturbations
class All_deformed : public Test_mesh
{
  cartdg::Solver sol;

  public:
  All_deformed()
  : sol{2, cartdg::config::max_row_size, 1.}
  {}

  virtual cartdg::Solver& solver() {return sol;}

  virtual std::vector<elem_handle> construct(cartdg::Boundary_condition* bc)
  {
    std::vector<elem_handle> handles;
    int bc_sn = sol.mesh().add_boundary_condition(bc);
    for (int i_elem = 0; i_elem < 9; ++i_elem) {
      std::vector<int> strides {3, 1};
      std::vector<int> inds;
      for (int i_dim = 0; i_dim < 2; ++i_dim) inds.push_back((i_elem/strides[i_dim])%3);
      int sn = sol.mesh().add_element(0, true, inds);
      handles.push_back({0, true, sn});
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        if (inds[i_dim] > 0) sol.mesh().connect_deformed(0, {handles[i_elem - strides[i_dim]].serial_n, sn}, {{i_dim, i_dim}, {1, 0}});
        if (inds[i_dim] != 1) sol.mesh().connect_boundary(0, true, sn, i_dim, (inds[i_dim] > 0), bc_sn);
      }
    }
    return handles;
  }
};

void test_marching(Test_mesh& tm, std::string name)
{
  // use `Copy` BCs. This is unstable for this case but it will still give the right answer as long as only one time step is executed
  auto handles = tm.construct(new cartdg::Copy);
  auto& sol = tm.solver();
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.initialize(Nonuniform_mass());
  sol.visualize_field(cartdg::State_variables(), "marching_" + name + "_init");
  // check that the iteration status is right at the start
  auto status = sol.iteration_status();
  REQUIRE(status.flow_time == 0.);
  REQUIRE(status.iteration == 0);
  // update
  sol.update();
  sol.visualize_field(cartdg::Physical_update(), "marching_" + name + "_diff");
  status = sol.iteration_status();
  REQUIRE(status.flow_time > 0.);
  REQUIRE(status.iteration == 1);
  // check that the computed update is approximately equal to the exact solution
  for (auto handle : handles) {
    for (int i_qpoint = 0; i_qpoint < sol.storage_params().n_qpoint(); ++i_qpoint) {
      for (int i_var = 0; i_var < sol.storage_params().n_var; ++i_var) {
        auto state   = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, cartdg::State_variables());
        auto update  = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, cartdg::Physical_update());
        auto correct = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, Nonuniform_residual());
        REQUIRE(update[i_var]/status.time_step == Approx(correct[i_var]).margin(1e-5*std::abs(state[i_var])));
      }
    }
  }
}

// test the solver on a sinusoid-derived initial condition which has a simple analytic solution
TEST_CASE("Solver time marching")
{
  SECTION("all cartesian")
  {
    All_cartesian ac;
    test_marching(ac, "car");
  }
  SECTION("all deformed")
  {
    All_deformed ad;
    test_marching(ad, "def");
  }
}
