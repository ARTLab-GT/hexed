#include <catch2/catch.hpp>
#include <config.hpp>
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
  virtual int n_var(int n_dim) const {return 4;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state) const
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0., 0.};
  }
};

class Normal_1 : public cartdg::Surface_func
{
  public:
  virtual int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const
  {
    return {outward_normal[1]};
  }
};

class Reciprocal_jacobian : public cartdg::Qpoint_func
{
  public:
  virtual int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(cartdg::Element& elem, const cartdg::Basis&, int i_qpoint, double time) const
  {
    return {1./elem.jacobian_determinant(i_qpoint)};
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

class Random_perturbation : public cartdg::Spacetime_func
{
  double rand_perturb() const
  {
    const int n = 1000;
    return 1. + 0.1/n*(rand()%n);
  };

  public:
  virtual int n_var(int n_dim) const {return n_dim + 2;};

  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n_dim {int(pos.size())};
    std::vector<double> result;
    double mass = 1.05;
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) result.push_back(mass*velocs[i_dim]*rand_perturb());
    result.push_back(mass*rand_perturb());
    result.push_back(1e5/0.4*rand_perturb());
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

class Parabola : public cartdg::Surface_geometry
{
  public:
  virtual std::array<double, 3> project_point(std::array<double, 3> point)
  {
    point[1] = 0.1*point[0]*point[0];
    return point;
  }
  virtual std::vector<double> line_intersections(std::array<double, 3> point0, std::array<double, 3> point1)
  {
    // only works for vertical lines
    if (std::abs(point0[0] - point1[0]) > 1e12) throw std::runtime_error("only for toy cases where the line is vertical");
    double intersection = 0.1*point0[0]*point0[0]/(point1[1] - point0[1]);
    return {10., intersection, -3.}; // add some other points just to try to confuse the face snapper
  }
};

TEST_CASE("Solver")
{
  static_assert (cartdg::config::max_row_size >= 3); // this test was written for row size 3
  cartdg::Solver sol {2, 3, 0.8};

  SECTION("initialization, sampling, and bounds")
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
    auto bounds = sol.bounds_field(cartdg::State_variables());
    REQUIRE(bounds.size() == 4);
    REQUIRE(bounds[0][0] == Approx(-.2*.2).scale(1.));
    REQUIRE(bounds[0][1] == Approx(1.6*.8).scale(1.));
    REQUIRE(bounds[2][0] == Approx(2.).scale(1.));
    REQUIRE(bounds[2][1] == Approx(2.).scale(1.));
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
    sol.mesh().connect_hanging(0, sn2, {sn3, sn4}, {{1, 1}, {0, 1}});
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

  SECTION("integrals")
  {
    SECTION("simple function, complex mesh")
    {
      int car0 = sol.mesh().add_element(0, false, {0, 0, 0});
      int car1 = sol.mesh().add_element(0, false, {1, 0, 0});
      int car2 = sol.mesh().add_element(1, false, {2, 2, 0});
      int sn0 = sol.mesh().add_element(0, true, {2, 0, 0});
      int sn1 = sol.mesh().add_element(0, true, {4, 0, 0});
      // connecting deformed elements with an empty space in between stretches them by a factor of 1.5
      sol.mesh().connect_deformed(0, {sn0, sn1}, {{0, 0}, {1, 0}});
      // add some boundary conditions for testing surface integrals
      int bc0 = sol.mesh().add_boundary_condition(new cartdg::Copy, new cartdg::Null_mbc);
      int bc1 = sol.mesh().add_boundary_condition(new cartdg::Copy, new cartdg::Null_mbc);
      int i = 0;
      for (int sn : {car0, car1, sn0, sn1}) {
        sol.mesh().connect_boundary(0, i++ >= 2, sn, 1, 0, bc0);
      }
      sol.mesh().connect_boundary(1, false, car2, 0, 1, bc0);
      sol.mesh().connect_boundary(1, false, car2, 1, 1, bc1);
      // finish setup
      sol.calc_jacobian();
      std::vector<double> state {0.3, -10., 0.7, 32.};
      sol.initialize(cartdg::Constant_func(state));
      // let's do some integrals
      auto integral = sol.integral_field(cartdg::State_variables());
      REQUIRE(integral.size() == 4);
      double area = 5.25*0.8*0.8;
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Approx(state[i_var]*area));
      }
      area = 5.5*0.8;
      integral = sol.integral_surface(cartdg::State_variables(), bc0);
      REQUIRE(integral.size() == 4);
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Approx(state[i_var]*area));
      }
      integral = sol.integral_surface(Normal_1(), bc0);
      REQUIRE(integral.size() == 1);
      REQUIRE(integral[0] == Approx(5*0.8));
      integral = sol.integral_surface(Normal_1(), bc1);
      REQUIRE(integral.size() == 1);
      REQUIRE(integral[0] == Approx(-0.5*0.8));
    }
    SECTION("complex function, simple mesh")
    {
      int sn = sol.mesh().add_element(0, false, {0, 0, 0});
      int bc0 = sol.mesh().add_boundary_condition(new cartdg::Nonpenetration, new cartdg::Null_mbc);
      sol.mesh().connect_boundary(0, false, sn, 0, 1, bc0);
      sol.calc_jacobian();
      sol.initialize(cartdg::Constant_func({0.3, 0., 0., 0.}));
      auto integral = sol.integral_field(Arbitrary_integrand());
      REQUIRE(integral.size() == 4);
      REQUIRE(integral[0] == Approx(std::pow(0.8, 3)/3*std::pow(0.8, 4)/4 - 0.8*0.8*0.3));
      integral = sol.integral_surface(Arbitrary_integrand(), bc0);
      REQUIRE(integral.size() == 4);
      REQUIRE(integral[0] == Approx(0.8*0.8*std::pow(0.8, 4)/4 - 0.8*0.3));
    }
  }

  SECTION("face/vertex snapping")
  {
    SECTION("Nominal_pos")
    {
      int el_sn = sol.mesh().add_element(1, true, {1, 2});
      int bc_sn = sol.mesh().add_boundary_condition(new cartdg::Copy, new cartdg::Nominal_pos);
      sol.mesh().connect_boundary(1, true, el_sn, 1, 1, bc_sn);
      sol.relax_vertices();
      sol.snap_vertices();
      sol.snap_faces();
      sol.calc_jacobian();
      // element should now be [(1 + 0.25)*0.8/2, (1 + 0.75)*0.8/2] x [(2 + 0.25)*0.8/2, 3*0.8/2]
      REQUIRE(sol.integral_field(cartdg::Constant_func({1.}))[0] == Approx(0.5*0.75*(0.8/2)*(0.8/2)));
    }
    SECTION("Surface_mbc")
    {
      int el_sn = sol.mesh().add_element(1, true, {0, 0});
      int bc_sn = sol.mesh().add_boundary_condition(new cartdg::Copy, new cartdg::Surface_mbc{new Parabola});
      sol.mesh().connect_boundary(1, true, el_sn, 1, 1, bc_sn);
      sol.snap_vertices();
      sol.calc_jacobian();
      // element should be a triangle
      REQUIRE(sol.integral_field(cartdg::Constant_func({1.}))[0] == Approx(0.5*(0.1*0.4*0.4)*0.4));
      sol.snap_faces();
      sol.calc_jacobian();
      // top element face should now be a parabola
      REQUIRE(sol.integral_field(cartdg::Constant_func({1.}))[0] == Approx(0.1*0.4*0.4*0.4/3.));
      // make sure that snapping again doesn't change anything
      sol.snap_faces();
      sol.calc_jacobian();
      REQUIRE(sol.integral_field(cartdg::Constant_func({1.}))[0] == Approx(0.1*0.4*0.4*0.4/3.));
    }
  }
}

class Test_mesh
{
  public:
  virtual cartdg::Solver& solver() = 0;
  virtual int bc_serial_n() = 0;
  virtual void construct(cartdg::Flow_bc* flow_bc) = 0;
};

// creates a 2x2x2 mesh
class All_cartesian : public Test_mesh
{
  int bc_sn;
  cartdg::Solver sol;
  public:

  All_cartesian()
  : sol{3, cartdg::config::max_row_size, 1.}
  {}

  virtual cartdg::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(cartdg::Flow_bc* flow_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, new cartdg::Null_mbc);
    std::vector<int> serial_n;
    for (int i_elem = 0; i_elem < 8; ++i_elem) {
      std::vector<int> strides {4, 2, 1};
      std::vector<int> inds;
      for (int i_dim = 0; i_dim < 3; ++i_dim) inds.push_back((i_elem/strides[i_dim])%2);
      int sn = sol.mesh().add_element(0, false, inds);
      serial_n.push_back(sn);
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        if (inds[i_dim]) sol.mesh().connect_cartesian(0, {serial_n[i_elem - strides[i_dim]], sn}, {i_dim});
        sol.mesh().connect_boundary(0, false, sn, i_dim, inds[i_dim], bc_sn);
      }
    }
  }
};

// creates a 3x3 mesh with the middle element rotated
class All_deformed : public Test_mesh
{
  cartdg::Solver sol;
  int bc_sn;
  bool rot_dir;

  public:
  All_deformed(bool rotation_direction)
  : sol{2, cartdg::config::max_row_size, 1.}, rot_dir{rotation_direction}
  {}

  virtual cartdg::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(cartdg::Flow_bc* flow_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, new cartdg::Null_mbc);
    std::vector<int> serial_n;
    for (int i_elem = 0; i_elem < 9; ++i_elem) {
      std::vector<int> strides {3, 1};
      std::vector<int> inds;
      for (int i_dim = 0; i_dim < 2; ++i_dim) inds.push_back((i_elem/strides[i_dim])%3);
      int sn = sol.mesh().add_element(0, true, inds);
      serial_n.push_back(sn);
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        if (inds[i_dim] > 0) {
          std::array<int, 2> i_elems {i_elem - strides[i_dim], i_elem};
          cartdg::Con_dir<cartdg::Deformed_element> dir {{i_dim, i_dim}, {1, 0}};
          for (int i_side = 0; i_side < 2; ++i_side) {
            if (i_elems[i_side] == 4) {
              dir.i_dim[i_side] = 1 - dir.i_dim[i_side];
              if (i_dim == rot_dir) dir.face_sign[i_side] = !dir.face_sign[i_side];
            }
          }
          sol.mesh().connect_deformed(0, {serial_n[i_elems[0]], serial_n[i_elems[1]]}, dir);
        }
        if (inds[i_dim] != 1) sol.mesh().connect_boundary(0, true, sn, i_dim, (inds[i_dim] > 0), bc_sn);
      }
    }
  }
};

// creates a mesh involving hanging vertices
class Hanging : public Test_mesh
{
  cartdg::Solver sol;
  int bc_sn;
  public:
  Hanging() : sol{2, cartdg::config::max_row_size, 1.} {}
  virtual cartdg::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}
  virtual void construct(cartdg::Flow_bc* flow_bc)
  {
    std::vector<int> serial_n;
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, new cartdg::Null_mbc);
    for (int i = 0; i < 2; ++i) {
      serial_n.push_back(sol.mesh().add_element(0, false, {i  , 0}));
      serial_n.push_back(sol.mesh().add_element(1, false, {i  , 2}));
      serial_n.push_back(sol.mesh().add_element(1,  true, {i+2, 2}));
    }
    sol.mesh().connect_cartesian(0, {serial_n[0], serial_n[3]}, {0});
    sol.mesh().connect_cartesian(1, {serial_n[1], serial_n[4]}, {0});
    sol.mesh().connect_deformed (1, {serial_n[2], serial_n[5]}, {{0, 0}, {1, 0}});
    sol.mesh().connect_cartesian(1, {serial_n[4], serial_n[2]}, {0}, {false, true});
    for (int i = 0; i < 2; ++i) {
      sol.mesh().connect_hanging(0, serial_n[3*i], {serial_n[1 + i], serial_n[4 + i]}, {{1, 1}, {1, 0}}, false, {bool(i), bool(i)});
      sol.mesh().connect_boundary(0, false, serial_n[3*i        ], 1, 0, bc_sn);
      sol.mesh().connect_boundary(0, false, serial_n[3*i        ], 0, i, bc_sn);
      sol.mesh().connect_boundary(1, false, serial_n[3*i + 1    ], 1, 1, bc_sn);
      sol.mesh().connect_boundary(1,  true, serial_n[3*i + 2    ], 1, 1, bc_sn);
      sol.mesh().connect_boundary(1,     i, serial_n[3*i + 1 + i], 0, i, bc_sn);
    }
    for (int i = 0; i < 2; ++i) sol.relax_vertices();
  }
};

// creates a single 3d element and then extrudes all faces
class Extrude_3d : public Test_mesh
{
  cartdg::Solver sol;
  int bc_sn;

  public:
  Extrude_3d()
  : sol{3, cartdg::config::max_row_size, .8}
  {}

  virtual cartdg::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(cartdg::Flow_bc* flow_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, new cartdg::Null_mbc);
    std::vector<cartdg::Mesh::elem_handle> handles;
    sol.mesh().add_element(0, true, {0, 0, 0});
    sol.mesh().extrude();
    sol.mesh().connect_rest(bc_sn);
  }
};

class Extrude_hanging : public Test_mesh
{
  cartdg::Solver sol;
  int bc_sn;

  public:
  Extrude_hanging()
  : sol{3, cartdg::config::max_row_size, .8}
  {}

  virtual cartdg::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(cartdg::Flow_bc* flow_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, new cartdg::Null_mbc);
    std::vector<cartdg::Mesh::elem_handle> handles;
    int coarse = sol.mesh().add_element(0, true, {-1, 0, 0});
    std::vector<int> fine(4);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        fine[2*i + j] = sol.mesh().add_element(1, true, {0, i, j});
        if (i) sol.mesh().connect_deformed(1, {fine[  j], fine[2   + j]}, {{1, 1}, {1, 0}});
        if (j) sol.mesh().connect_deformed(1, {fine[2*i], fine[2*i + 1]}, {{2, 2}, {1, 0}});
      }
    }
    int extra [2];
    for (int i = 0; i < 2; ++i) {
      extra[i] = sol.mesh().add_element(1, true, {0, 2, i});
      sol.mesh().connect_deformed(1, {fine[2 + i], extra[i]}, {{1, 1}, {1, 0}});
    }
    sol.mesh().connect_deformed(1, {extra[0], extra[1]}, {{2, 2}, {1, 0}});
    sol.mesh().connect_hanging(0, coarse, fine, {{0, 0}, {1, 0}}, true, {true, true, true, true});
    sol.mesh().extrude();
    sol.mesh().connect_rest(bc_sn);
    for (int i = 0; i < 2; ++i) sol.relax_vertices();
  }
};

void test_marching(Test_mesh& tm, std::string name)
{
  // use `Copy` BCs. This is unstable for this case but it will still give the right answer as long as only one time step is executed
  tm.construct(new cartdg::Copy);
  auto& sol = tm.solver();
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.initialize(Nonuniform_mass());
  #if CARTDG_USE_TECPLOT
  sol.visualize_field_tecplot(cartdg::State_variables(), "marching_" + name + "_init");
  sol.visualize_surface_tecplot(tm.bc_serial_n(), "marching_" + name + "_surface");
  #endif
  // check that the iteration status is right at the start
  auto status = sol.iteration_status();
  REQUIRE(status.flow_time == 0.);
  REQUIRE(status.iteration == 0);
  // update
  sol.update();
  #if CARTDG_USE_TECPLOT
  sol.visualize_field_tecplot(cartdg::Physical_update(), "marching_" + name + "_diff");
  #endif
  status = sol.iteration_status();
  REQUIRE(status.flow_time > 0.);
  REQUIRE(status.iteration == 1);
  // check that the computed update is approximately equal to the exact solution
  for (auto handle : sol.mesh().elem_handles()) {
    for (int i_qpoint = 0; i_qpoint < sol.storage_params().n_qpoint(); ++i_qpoint) {
      for (int i_var = 0; i_var < sol.storage_params().n_var; ++i_var) {
        auto state   = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, cartdg::State_variables());
        auto update  = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, cartdg::Physical_update());
        auto correct = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, Nonuniform_residual());
        REQUIRE(update[i_var]/status.time_step == Approx(correct[i_var]).margin(1e-3*std::abs(state[i_var])));
      }
    }
  }
}

void test_conservation(Test_mesh& tm, std::string name)
{
  srand(406);
  tm.construct(new cartdg::Nonpenetration());
  auto& sol = tm.solver();
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.initialize(Random_perturbation());
  // check that the iteration status is right at the start
  auto status = sol.iteration_status();
  // update
  sol.update();
  status = sol.iteration_status();
  // check that the computed update is approximately equal to the exact solution
  for (int i_var : {sol.storage_params().n_var - 2, sol.storage_params().n_var - 1}) {
    auto state  = sol.integral_field(cartdg::State_variables());
    auto update = sol.integral_field(cartdg::Physical_update());
    REQUIRE(update[i_var]/status.time_step == Approx(0.).scale(std::abs(state[i_var])));
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
    All_deformed ad0 (0);
    test_marching(ad0, "def0");
    All_deformed ad1 (1);
    test_marching(ad1, "def1");
  }
  SECTION("with hanging nodes")
  {
    Hanging hg;
    test_marching(hg, "hanging");
  }
  SECTION("3d extruded")
  {
    Extrude_3d e3;
    test_marching(e3, "extrude_3d");
  }
  SECTION("extruded with deformed hanging nodes")
  {
    Extrude_hanging eh;
    #if CARTDG_USE_OTTER
    SECTION("vis")
    {
      eh.construct(new cartdg::Copy);
      eh.solver().calc_jacobian();
      otter::plot plt;
      eh.solver().visualize_edges_otter(plt);
      eh.solver().visualize_surface_otter(plt, eh.bc_serial_n(), otter::const_colormap(otter::colors::css4["darkgrey"]), cartdg::Pressure(), {0., 1.}, true);
      plt.show();
      eh.solver().mesh().valid().assert_valid();
    }
    #endif
    //test_marching(eh, "extrude_hanging");
  }
}

// test the solver on a randomly perturbed input (for which it can't possibly be accurate) and verify conservation
TEST_CASE("Solver conservation")
{
  SECTION("all cartesian")
  {
    All_cartesian ac;
    test_conservation(ac, "car");
  }
  SECTION("all deformed")
  {
    All_deformed ad0 (0);
    test_conservation(ad0, "def0");
    All_deformed ad1 (1);
    test_conservation(ad1, "def1");
  }
  SECTION("with hanging nodes")
  {
    Hanging hg;
    test_conservation(hg, "hanging");
  }
  SECTION("3d extruded")
  {
    Extrude_3d e3;
    test_conservation(e3, "extrude_3d");
  }
  SECTION("extruded with deformed hanging nodes")
  {
    Extrude_hanging eh;
    test_conservation(eh, "extrude_hanging");
  }
}

TEST_CASE("face extrusion")
{
  SECTION("2D")
  {
    int serial_n [3][3];
    serial_n[1][1] = -1; // so that we know if we accidentally use this
    cartdg::Solver solver {2, 2, 1.};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if ((i != 1) || (j != 1)) {
          serial_n[i][j] = solver.mesh().add_element(1, true, {i, j});
          if ((i > 0) && (j != 1)) solver.mesh().connect_deformed(1, {serial_n[i-1][j], serial_n[i][j]}, {{0, 0}, {1, 0}});
          if ((j > 0) && (i != 1)) solver.mesh().connect_deformed(1, {serial_n[i][j-1], serial_n[i][j]}, {{1, 1}, {1, 0}});
        }
      }
    }
    solver.mesh().extrude();
    solver.calc_jacobian();
    REQUIRE(solver.integral_field(Reciprocal_jacobian())[0] == Approx(24./4.)); // check number of elements
    for (int i = 0; i < 3; ++i) solver.relax_vertices(); // so that we can see better
    auto valid = solver.mesh().valid();
    REQUIRE(valid.n_duplicate == 0);
    REQUIRE(valid.n_missing == 16);
  }
  SECTION("3D")
  {
    int serial_n [3][3][3];
    serial_n[1][1][1] = -1; // so that we know if we accidentally use this
    cartdg::Solver solver {3, 2, 1.};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          if ((i != 1) || (j != 1) || (k != 1)) {
            serial_n[i][j][k] = solver.mesh().add_element(0, true, {i, j, k});
            if ((i > 0) && ((j != 1) || (k != 1))) solver.mesh().connect_deformed(0, {serial_n[i-1][j][k], serial_n[i][j][k]}, {{0, 0}, {1, 0}});
            if ((j > 0) && ((i != 1) || (k != 1))) solver.mesh().connect_deformed(0, {serial_n[i][j-1][k], serial_n[i][j][k]}, {{1, 1}, {1, 0}});
            if ((k > 0) && ((i != 1) || (j != 1))) solver.mesh().connect_deformed(0, {serial_n[i][j][k-1], serial_n[i][j][k]}, {{2, 2}, {1, 0}});
          }
        }
      }
    }
    solver.mesh().extrude();
    solver.calc_jacobian();
    REQUIRE(solver.integral_field(Reciprocal_jacobian())[0] == Approx(86.)); // check number of elements
    for (int i = 0; i < 3; ++i) solver.relax_vertices(); // so that we can see better
    auto valid = solver.mesh().valid();
    REQUIRE(valid.n_duplicate == 0);
    REQUIRE(valid.n_missing == 60);
    #if CARTDG_USE_TECPLOT
    solver.visualize_field_tecplot(cartdg::State_variables(), "extrude");
    #endif
  }
}
