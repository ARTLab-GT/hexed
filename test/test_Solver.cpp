#include <catch2/catch_all.hpp>
#include <hexed/config.hpp>
#include <hexed/Solver.hpp>

class Arbitrary_initializer : public hexed::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 4;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {pos[0]*pos[1], 1., 2., 3.};
  }
};

class Bad_initializer : public hexed::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    return {0., 0.};
  }
};

class Arbitrary_integrand : public hexed::Domain_func
{
  public:
  virtual int n_var(int n_dim) const {return 4;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state) const
  {
    return std::vector<double> {pos[0]*pos[0]*pos[1]*pos[1]*pos[1] - state[0] + time, 0., 0., 0.};
  }
};

class Normal_1 : public hexed::Surface_func
{
  public:
  virtual int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time,
                                         std::vector<double> state, std::vector<double> outward_normal) const
  {
    return {outward_normal[1]};
  }
};

class Reciprocal_jacobian : public hexed::Qpoint_func
{
  public:
  virtual int n_var(int n_dim) const {return 1;}
  virtual std::vector<double> operator()(hexed::Element& elem, const hexed::Basis&, int i_qpoint, double time) const
  {
    return {1./elem.jacobian_determinant(i_qpoint)};
  }
};

class Sinusoid_veloc0 : public hexed::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return n_dim + 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    // constant density and pressure, but sinusoidally-varying velocity with zero divergence
    const int n_dim = pos.size();
    double cos_arg = 0;
    for (double p : pos) cos_arg += p;
    double cos_part = std::cos(cos_arg);
    std::vector<double> state(n_var(n_dim), 0.);
    for (int i_dim = 0; i_dim < n_dim; ++i_dim) state[i_dim] = 1.2*cos_part;
    state[0] *= 1 - n_dim;
    state[n_dim] = 1.2;
    state[n_dim + 1] = 1e5/0.4 + .5*1.2*(n_dim - 1 + (1 - n_dim)*(1 - n_dim))*cos_part*cos_part;
    return state;
  }
};

class Sinusoid_veloc1 : public hexed::Spacetime_func
{
  public:
  virtual int n_var(int n_dim) const {return n_dim + 2;}
  // requires `n_dim >= 2`
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n_dim = pos.size();
    std::vector<double> state(n_var(n_dim), 0.);
    state[n_dim] = 2.3;
    state[0] =  state[n_dim]*.1*std::sin(pos[0]);
    state[1] = -state[n_dim]*.2*std::cos(pos[1]);
    state[n_dim + 1] = 10.;
    return state;
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

class Nonuniform_mass : public hexed::Spacetime_func
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

class Random_perturbation : public hexed::Spacetime_func
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

class Nonuniform_residual : public hexed::Spacetime_func
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

class Tanh : public hexed::Spacetime_func
{
  double scale;
  public:
  Tanh(double scale_arg) : scale{scale_arg} {}
  virtual int n_var(int n_dim) const {return n_dim + 2;}
  virtual std::vector<double> operator()(std::vector<double> pos, double time) const
  {
    const int n_dim {int(pos.size())};
    std::vector<double> result(n_dim + 2, 0.);
    result[0] = 1. - .1*std::tanh((pos[0] - .5)/scale);
    result[n_dim] = 1.;
    result[n_dim + 1] = 1.5;
    return result;
  }
};

class Parabola : public hexed::Surface_geom
{
  public:
  hexed::Nearest_point<Eigen::Dynamic> nearest_point(hexed::Mat<> point, double, double) override
  {
    return hexed::Nearest_point<Eigen::Dynamic>(point); // empty point
  }
  std::vector<double> intersections(hexed::Mat<> point0, hexed::Mat<> point1) override
  {
    // only works for vertical lines
    if (std::abs(point0[0] - point1[0]) > 1e12) throw std::runtime_error("only for toy cases where the line is vertical");
    double intersection = 1. + 0.3*point0[0]*point0[0]/(point1[1] - point0[1]);
    return {10., intersection, -3.}; // add some other points just to try to confuse the face snapper
  }
};

class Shrink_pos0 : public hexed::Mesh_bc
{
  public:
  virtual void snap_vertices(hexed::Boundary_connection& con)
  {
    const int stride = hexed::math::pow(2, con.storage_params().n_dim - 1 - con.i_dim());
    for (int i_vert = 0; i_vert < con.storage_params().n_vertices(); ++i_vert) {
      if ((i_vert/stride)%2 == con.inside_face_sign()) {
        con.element().vertex(i_vert).pos[0] = .1*.8;
      }
    }
  }
  virtual void snap_node_adj(hexed::Boundary_connection&, const hexed::Basis&) {}
};

class Boundary_perturbation : public hexed::Mesh_bc
{
  double arbitrary [12] {.5, -.2, -.7, .6, .8, -.8, .1, .3, -.9, -.1, .4, .6};
  public:
  virtual void snap_vertices(hexed::Boundary_connection&) {}
  virtual void snap_node_adj(hexed::Boundary_connection& con, const hexed::Basis&)
  {
    if (!con.element().node_adjustments()) return; // Cartesian elements don't have `node_adjustments()`, so in this case just exit
    auto params {con.storage_params()};
    const int nfq = params.n_qpoint()/params.row_size;
    for (int i_qpoint = 0; i_qpoint < nfq; ++i_qpoint) {
      con.element().node_adjustments()[(2*con.i_dim() + con.inside_face_sign())*nfq + i_qpoint] += 0.02*arbitrary[i_qpoint%12];
    }
  }
};

TEST_CASE("Solver")
{
  static_assert (hexed::config::max_row_size >= 3); // this test was written for row size 3
  hexed::Solver sol {2, 3, 0.8, true};
  int catchall_bc = sol.mesh().add_boundary_condition(new hexed::Copy(), new hexed::Null_mbc());

  SECTION("initialization and introspection")
  {
    sol.mesh().add_element(0, false, {0, 0, 0});
    int sn0 = sol.mesh().add_element(0, false, {1, 0, 0});
    int sn1 = sol.mesh().add_element(2, true, {-1, 0, 0});
    sol.mesh().cleanup();
    // initialization/sampling
    REQUIRE_THROWS(sol.initialize(Bad_initializer())); // if number of variables of func is wrong, should throw
    REQUIRE_THROWS(sol.initialize(Arbitrary_initializer())); // mesh must be valid before you can initialize
    sol.mesh().connect_rest(catchall_bc);
    sol.mesh().cleanup();
    sol.initialize(Arbitrary_initializer());
    auto sample = sol.sample(0, false, sn0, 4, hexed::State_variables()); // sample the midpoint of the element because we know the exact position
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Catch::Approx(1.2*0.4));
    REQUIRE(sample[1] == Catch::Approx(1.));
    REQUIRE(sample[2] == Catch::Approx(2.));
    REQUIRE(sample[3] == Catch::Approx(3.));
    sample = sol.sample(2, true, sn1, 4, hexed::State_variables());
    REQUIRE(sample.size() == 4);
    REQUIRE(sample[0] == Catch::Approx(-0.1*0.1));
    // variable bounds calculation
    auto bounds = sol.bounds_field(hexed::State_variables());
    REQUIRE(bounds.size() == 4);
    REQUIRE(bounds[0][0] == Catch::Approx(-.2*.2).scale(1.));
    REQUIRE(bounds[0][1] == Catch::Approx(1.6*.8).scale(1.));
    REQUIRE(bounds[2][0] == Catch::Approx(2.).scale(1.));
    REQUIRE(bounds[2][1] == Catch::Approx(2.).scale(1.));
    // artificial viscosity setting
    sol.set_art_visc_constant(.28);
    REQUIRE(sol.sample(0, false, sn0, 4, hexed::Art_visc_coef())[0] == Catch::Approx(.28));
  }

  SECTION("vertex relaxation")
  {
    int sn0 = sol.mesh().add_element(0, false, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0,  true, {1, 0, 0});
    sol.mesh().connect_cartesian(0, {sn0, sn1}, {0}, {false, true});
    sol.mesh().cleanup();
    sol.mesh().relax(0.5);
    REQUIRE(sol.sample(0, false, sn0, 4, hexed::Position_func())[0] == Catch::Approx(0.8*0.5));
    REQUIRE(sol.sample(0,  true, sn1, 4, hexed::Position_func())[0] == Catch::Approx(0.8*1.375));
    REQUIRE(sol.sample(0,  true, sn1, 4, hexed::Position_func())[1] == Catch::Approx(0.8*0.5));
  }

  SECTION("local time step scale")
  {
    int sn0 = sol.mesh().add_element(0,  true, {0, 0, 0});
    int sn1 = sol.mesh().add_element(0, false, {-1, 0, 0});
    int sn2 = sol.mesh().add_element(1, false, {-2, -1, 0});
    int sn3 = sol.mesh().add_element(1, false, {-1, -1, 0});
    sol.mesh().connect_cartesian(0, {sn1, sn0}, {0}, {false, true});
    sol.mesh().connect_hanging(0, sn1, {sn2, sn3}, {{1, 1}, {0, 1}});
    int bc = sol.mesh().add_boundary_condition(new hexed::Copy, new Shrink_pos0);
    sol.mesh().connect_boundary(0, true, sn0, 0, 1, bc);
    sol.mesh().connect_rest(catchall_bc);
    sol.mesh().cleanup();
    sol.mesh().valid().assert_valid();
    sol.mesh().relax(0.);
    sol.calc_jacobian();
    sol.initialize(hexed::Constant_func({0., 0., 1.4, 1./.4})); // initialize with zero velocity and unit speed of sound
    sol.update(); // sets the time step scale (among other things)
    double cfl = hexed::Gauss_legendre(3).max_cfl()/2.; // divide by 2 cause that's the number of dimensions
    // in sn0, TSS is max_cfl*root_size*determinant/normal_sum/sound_speed
    REQUIRE(sol.sample(0,  true, sn0, 4, hexed::Time_step_scale_func())[0] == Catch::Approx(0.7*cfl*.8*2/11./1.));
    // in sn1, TSS varies bilinearly
    REQUIRE(sol.sample(0, false, sn1, 4, hexed::Time_step_scale_func())[0] == Catch::Approx(0.7*cfl*.8*(2*2/11. + .5 + 1.)/4.));
    // TSS at hanging nodes should be set to match coarse element
    REQUIRE(sol.sample(1, false, sn2, 4, hexed::Time_step_scale_func())[0] == Catch::Approx(0.7*cfl*.8*(3*.5 + (.5 + 2/11.)/2.)/4.));
    REQUIRE(sol.sample(1, false, sn3, 4, hexed::Time_step_scale_func())[0] == Catch::Approx(0.7*cfl*.8*(2*.5 + (.5 + 2/11.)/2. + 2/11.)/4.));
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
      int bc0 = sol.mesh().add_boundary_condition(new hexed::Copy, new hexed::Null_mbc);
      int bc1 = sol.mesh().add_boundary_condition(new hexed::Copy, new hexed::Null_mbc);
      int bc2 = sol.mesh().add_boundary_condition(new hexed::Copy, new hexed::Null_mbc);
      int i = 0;
      for (int sn : {car0, car1, sn0, sn1}) {
        sol.mesh().connect_boundary(0, i++ >= 2, sn, 1, 0, bc0);
      }
      sol.mesh().connect_boundary(1, false, car2, 0, 1, bc0);
      sol.mesh().connect_boundary(1, false, car2, 1, 1, bc1);
      sol.mesh().connect_rest(bc2);
      // finish setup
      sol.mesh().cleanup();
      sol.calc_jacobian();
      std::vector<double> state {0.3, -10., 0.7, 32.};
      sol.initialize(hexed::Constant_func(state));
      // let's do some integrals
      auto integral = sol.integral_field(hexed::State_variables());
      REQUIRE(integral.size() == 4);
      double area = 5.25*0.8*0.8;
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Catch::Approx(state[i_var]*area));
      }
      area = 5.5*0.8;
      integral = sol.integral_surface(hexed::State_variables(), bc0);
      REQUIRE(integral.size() == 4);
      for (int i_var = 0; i_var < 4; ++i_var) {
        REQUIRE(integral[i_var] == Catch::Approx(state[i_var]*area));
      }
      integral = sol.integral_surface(Normal_1(), bc0);
      REQUIRE(integral.size() == 1);
      REQUIRE(integral[0] == Catch::Approx(5*0.8));
      integral = sol.integral_surface(Normal_1(), bc1);
      REQUIRE(integral.size() == 1);
      REQUIRE(integral[0] == Catch::Approx(-0.5*0.8));
    }
    SECTION("complex function, simple mesh")
    {
      int sn = sol.mesh().add_element(0, false, {0, 0, 0});
      int bc0 = sol.mesh().add_boundary_condition(new hexed::Nonpenetration, new hexed::Null_mbc);
      int bc1 = sol.mesh().add_boundary_condition(new hexed::Copy, new hexed::Null_mbc);
      sol.mesh().connect_boundary(0, false, sn, 0, 1, bc0);
      sol.mesh().connect_rest(bc1);
      sol.mesh().cleanup();
      sol.calc_jacobian();
      sol.initialize(hexed::Constant_func({0.3, 0., 0., 0.}));
      auto integral = sol.integral_field(Arbitrary_integrand());
      REQUIRE(integral.size() == 4);
      REQUIRE(integral[0] == Catch::Approx(std::pow(0.8, 3)/3*std::pow(0.8, 4)/4 - 0.8*0.8*0.3));
      integral = sol.integral_surface(Arbitrary_integrand(), bc0);
      REQUIRE(integral.size() == 4);
      REQUIRE(integral[0] == Catch::Approx(0.8*0.8*std::pow(0.8, 4)/4 - 0.8*0.3));
    }
  }

  SECTION("face/vertex snapping")
  {
    SECTION("Nominal_pos")
    {
      int el_sn = sol.mesh().add_element(1, true, {1, 2});
      int bc_sn = sol.mesh().add_boundary_condition(new hexed::Copy, new hexed::Nominal_pos);
      sol.mesh().connect_boundary(1, true, el_sn, 1, 1, bc_sn);
      sol.mesh().connect_rest(catchall_bc);
      sol.mesh().cleanup();
      sol.mesh().relax(.5);
      sol.mesh().valid().assert_valid();
      sol.calc_jacobian();
      // element should now be [(1 + 0.25)*0.8/2, (1 + 0.75)*0.8/2] x [(2 + 0.25)*0.8/2, 3*0.8/2]
      REQUIRE(sol.integral_field(hexed::Constant_func({1.}))[0] == Catch::Approx(0.5*0.75*(0.8/2)*(0.8/2)));
    }
    SECTION("Surface_mbc")
    {
      static_assert (hexed::config::max_row_size >= 4);
      hexed::Solver sol4 {2, 4, 0.8, true};
      int catchall_bc4 = sol4.mesh().add_boundary_condition(new hexed::Copy(), new hexed::Null_mbc());
      int el_sn = sol4.mesh().add_element(1, true, {0, 0});
      int bc_sn = sol4.mesh().add_boundary_condition(new hexed::Copy, new hexed::Geom_mbc{new Parabola});
      sol4.mesh().connect_boundary(1, true, el_sn, 1, 1, bc_sn);
      sol4.mesh().connect_rest(catchall_bc4);
      sol4.mesh().cleanup();
      sol4.mesh().valid().assert_valid();
      sol4.calc_jacobian();
      // top element face should now be a parabola
      REQUIRE(sol4.integral_field(hexed::Constant_func({1.}))[0] - .4*.4 == Catch::Approx(.3*0.4*0.4*0.4/3.));
      // make sure that snapping again doesn't change anything
      sol4.calc_jacobian();
      REQUIRE(sol4.integral_field(hexed::Constant_func({1.}))[0] - .4*.4 == Catch::Approx(.3*0.4*0.4*0.4/3.));
    }
  }
}

class Test_mesh
{
  public:
  virtual hexed::Solver& solver() = 0;
  virtual int bc_serial_n() = 0;
  virtual void construct(hexed::Flow_bc* flow_bc, hexed::Mesh_bc* mesh_bc) = 0;
};

// creates a 2x2x2 mesh
class All_cartesian : public Test_mesh
{
  int bc_sn;
  hexed::Solver sol;
  public:

  All_cartesian(bool local, hexed::Transport_model transport = hexed::inviscid)
  : sol{3, hexed::config::max_row_size, 1., local, transport}
  {}

  virtual hexed::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(hexed::Flow_bc* flow_bc, hexed::Mesh_bc* mesh_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, mesh_bc);
    std::vector<int> serial_n;
    for (int i_elem = 0; i_elem < 8; ++i_elem) {
      std::vector<int> strides {4, 2, 1};
      std::vector<int> inds;
      for (int i_dim = 0; i_dim < 3; ++i_dim) inds.push_back((i_elem/strides[i_dim])%2);
      int sn = sol.mesh().add_element(0, false, inds);
      serial_n.push_back(sn);
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        if (inds[i_dim]) sol.mesh().connect_cartesian(0, {serial_n[i_elem - strides[i_dim]], sn}, {i_dim});
        //if (inds[i_dim]) sol.mesh().connect_cartesian(0, {sn, serial_n[i_elem - strides[i_dim]]}, {i_dim});
        sol.mesh().connect_boundary(0, false, sn, i_dim, inds[i_dim], bc_sn);
      }
    }
    sol.mesh().cleanup();
  }
};

// creates a 3x3 mesh with the middle element rotated
class All_deformed : public Test_mesh
{
  hexed::Solver sol;
  int bc_sn;
  bool rot_dir;

  public:
  All_deformed(bool rotation_direction, bool local, hexed::Transport_model transport = hexed::inviscid)
  : sol{2, hexed::config::max_row_size, 1., local, transport}, rot_dir{rotation_direction}
  {}

  virtual hexed::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(hexed::Flow_bc* flow_bc, hexed::Mesh_bc* mesh_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, mesh_bc);
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
          hexed::Con_dir<hexed::Deformed_element> dir {{i_dim, i_dim}, {1, 0}};
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
    sol.mesh().cleanup();
  }
};

class Extrude_hanging : public Test_mesh
{
  hexed::Solver sol;
  int bc_sn;
  int id;
  int jd;
  int kd;

  public:
  Extrude_hanging(int i_dim, int j_dim, bool local, hexed::Transport_model transport = hexed::inviscid)
  : sol{3, hexed::config::max_row_size, .8, local, transport}, id{i_dim}, jd{j_dim}, kd{3 - i_dim - j_dim}
  {}

  virtual hexed::Solver& solver() {return sol;}
  virtual int bc_serial_n() {return bc_sn;}

  virtual void construct(hexed::Flow_bc* flow_bc, hexed::Mesh_bc* mesh_bc)
  {
    bc_sn = sol.mesh().add_boundary_condition(flow_bc, mesh_bc);
    std::vector<hexed::Mesh::elem_handle> handles;
    std::vector<int> coarse_pos(3, 0);
    coarse_pos[id] = -1;
    int coarse = sol.mesh().add_element(0, true, coarse_pos);
    std::vector<int> fine(4);
    int dim0 = std::min(jd, kd);
    int dim1 = std::max(jd, kd);
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        std::vector<int> fine_pos(3, 0);
        fine_pos[dim0] = i;
        fine_pos[dim1] = j;
        fine[2*i + j] = sol.mesh().add_element(1, true, fine_pos);
        if (i) sol.mesh().connect_deformed(1, {fine[  j], fine[2   + j]}, {{dim0, dim0}, {1, 0}});
        if (j) sol.mesh().connect_deformed(1, {fine[2*i], fine[2*i + 1]}, {{dim1, dim1}, {1, 0}});
      }
    }
    int extra [2];
    for (int i = 0; i < 2; ++i) {
      std::vector<int> fine_pos(3, 0);
      fine_pos[jd] = 2;
      fine_pos[kd] = i;
      extra[i] = sol.mesh().add_element(1, true, fine_pos);
      sol.mesh().connect_deformed(1, {fine[1 + (jd < kd) + i*(1 + (jd > kd))], extra[i]}, {{jd, jd}, {1, 0}});
    }
    sol.mesh().connect_deformed(1, {extra[0], extra[1]}, {{kd, kd}, {1, 0}});
    sol.mesh().connect_hanging(0, coarse, fine, {{id, id}, {1, 0}}, true, {true, true, true, true});
    sol.mesh().extrude();
    sol.mesh().connect_rest(bc_sn);
    sol.mesh().cleanup();
    for (int i = 0; i < 2; ++i) sol.mesh().relax(.5);
  }
};

void test_marching(Test_mesh& tm, std::string name)
{
  // use `Copy` BCs. This is unstable for this case but it will still give the right answer as long as only one time step is executed
  tm.construct(new hexed::Copy, new hexed::Null_mbc);
  auto& sol = tm.solver();
  sol.nspace().assign("use_filter", 0);
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.initialize(Nonuniform_mass());
  // check that the iteration status is right at the start
  auto status = sol.iteration_status();
  REQUIRE(status.flow_time == 0.);
  REQUIRE(status.iteration == 0);
  // update
  sol.nspace().assign("max_safety", 1e-3);
  sol.update();
  sol.compute_residual();
  status = sol.iteration_status();
  REQUIRE(status.flow_time > 0.);
  REQUIRE(status.iteration == 1);
  // check that the computed update is approximately equal to the exact solution
  for (auto handle : sol.mesh().elem_handles()) {
    for (int i_qpoint = 0; i_qpoint < sol.storage_params().n_qpoint(); ++i_qpoint) {
      for (int i_var = 0; i_var < sol.storage_params().n_var; ++i_var) {
        auto state   = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, hexed::State_variables());
        auto update  = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, hexed::Physical_residual());
        auto correct = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, Nonuniform_residual());
        REQUIRE(update[i_var]/status.time_step == Catch::Approx(correct[i_var]).margin(1e-3*std::abs(state[i_var])));
      }
    }
  }
}

void test_visc(Test_mesh& tm, std::string name)
{
  // use `Copy` BCs. This is unstable for this case but it will still give the right answer as long as only one time step is executed
  tm.construct(new hexed::Copy, new hexed::Null_mbc);
  auto& sol = tm.solver();
  const int n_dim = sol.storage_params().n_dim;
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.initialize(Sinusoid_veloc0());
  // update
  sol.nspace().assign("max_safety", 1e-4);
  sol.update();
  sol.compute_residual();
  auto status = sol.iteration_status();
  // check that the computed update is approximately equal to the exact solution
  for (auto handle : sol.mesh().elem_handles()) {
    for (int i_qpoint = 0; i_qpoint < sol.storage_params().n_qpoint(); ++i_qpoint) {
      auto state  = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, hexed::State_variables());
      auto update = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, hexed::Physical_residual());
      auto tss    = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, i_qpoint, hexed::Time_step_scale_func());
      for (int i_var = 0; i_var < n_dim; ++i_var) {
        double margin = 1.2*(hexed::config::max_row_size > 6 ? 1e-3 : 1.);
        // rate of change of momentum should be sinusoidal and rate of change of scalars should be 0
        REQUIRE(update[i_var]/status.time_step == Catch::Approx(-n_dim*3.*(state[i_var] - update[i_var]*tss[0])/1.2).margin(margin));
      }
    }
  }
}

void test_conservation(Test_mesh& tm, std::string name)
{
  srand(406);
  tm.construct(new hexed::Nonpenetration(), new Boundary_perturbation);
  auto& sol = tm.solver();
  sol.mesh().valid().assert_valid();
  sol.snap_faces();
  sol.calc_jacobian();
  sol.initialize(Random_perturbation());
  sol.set_art_visc_constant(30.);
  // check that the iteration status is right at the start
  auto status = sol.iteration_status();
  // update
  sol.nspace().assign("max_safety", .01);
  sol.update();
  sol.compute_residual();
  status = sol.iteration_status();
  auto state  = sol.integral_field(hexed::State_variables());
  auto update = sol.integral_field(hexed::Physical_residual());
  for (int i_var : {sol.storage_params().n_var - 2, sol.storage_params().n_var - 1}) {
    REQUIRE(update[i_var]/status.time_step == Catch::Approx(0.).scale(std::abs(state[i_var])));
  }
}

// test the solver on a sinusoid-derived initial condition which has a simple analytic solution
TEST_CASE("Solver time marching")
{
  SECTION("all cartesian")
  {
    All_cartesian ac(true);
    test_marching(ac, "car");
  }
  SECTION("all deformed")
  {
    All_deformed ad0 (0, true);
    test_marching(ad0, "def0");
    All_deformed ad1 (1, true);
    test_marching(ad1, "def1");
  }
  SECTION("extruded with deformed hanging nodes")
  {
    #define TEST_DIMENSIONS(i_dim, j_dim) \
      SECTION("dimensions " #i_dim " " #j_dim) { \
          Extrude_hanging eh(i_dim, j_dim, true); \
          test_marching(eh, "extrude_hanging"); \
      }
    #if NDEBUG
    TEST_DIMENSIONS(0, 1)
    TEST_DIMENSIONS(0, 2)
    TEST_DIMENSIONS(1, 0)
    TEST_DIMENSIONS(1, 2)
    TEST_DIMENSIONS(2, 0)
    TEST_DIMENSIONS(2, 1)
    #endif
    #undef TEST_DIMENSIONS
  }
}

// test the solver on a sinusoid-derived initial condition which has a simple analytic solution
TEST_CASE("Solver viscosity")
{
  SECTION("all cartesian")
  {
    All_cartesian ac(true, hexed::Transport_model::constant(3.));
    test_visc(ac, "car");
  }
  SECTION("all deformed")
  {
    All_deformed ad0 (0, true, hexed::Transport_model::constant(3.));
    test_visc(ad0, "def0");
    All_deformed ad1 (1, true, hexed::Transport_model::constant(3.));
    test_visc(ad1, "def1");
  }
  SECTION("extruded with deformed hanging nodes")
  {
    #define TEST_DIMENSIONS(i_dim, j_dim) \
      SECTION("dimensions " #i_dim " " #j_dim) { \
          Extrude_hanging eh(i_dim, j_dim, true, hexed::Transport_model::constant(3.)); \
          test_visc(eh, "extrude_hanging"); \
      }
    #if NDEBUG
    TEST_DIMENSIONS(0, 1)
    TEST_DIMENSIONS(0, 2)
    TEST_DIMENSIONS(1, 0)
    TEST_DIMENSIONS(1, 2)
    TEST_DIMENSIONS(2, 0)
    TEST_DIMENSIONS(2, 1)
    #endif
    #undef TEST_DIMENSIONS
  }
}

// test the solver on a randomly perturbed input (for which it can't possibly be accurate) and verify conservation
TEST_CASE("Solver conservation")
{
  SECTION("all cartesian")
  {
    All_cartesian ac(true);
    test_conservation(ac, "car");
  }
  SECTION("all deformed")
  {
    All_deformed ad0 (0, true);
    test_conservation(ad0, "def0");
    All_deformed ad1 (1, true);
    test_conservation(ad1, "def1");
  }
  SECTION("extruded with deformed hanging nodes")
  {
    #define TEST_CONSERVATION(i_dim, j_dim) \
      SECTION("dimensions " #i_dim " " #j_dim) { \
          Extrude_hanging eh(i_dim, j_dim, true); \
          test_conservation(eh, "extrude_hanging"); \
      } \

    #if NDEBUG
    TEST_CONSERVATION(0, 1)
    TEST_CONSERVATION(0, 2)
    TEST_CONSERVATION(1, 0)
    TEST_CONSERVATION(1, 2)
    TEST_CONSERVATION(2, 0)
    TEST_CONSERVATION(2, 1)
    #endif
  }
}

TEST_CASE("face extrusion")
{
  SECTION("2D")
  {
    int serial_n [3][3];
    serial_n[1][1] = -1; // so that we know if we accidentally use this
    hexed::Solver solver {2, 2, 1.};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if ((i != 1) || (j != 1)) {
          serial_n[i][j] = solver.mesh().add_element(1, true, {i, j});
          if ((i > 0) && (j != 1)) solver.mesh().connect_deformed(1, {serial_n[i-1][j], serial_n[i][j]}, {{0, 0}, {1, 0}});
          if ((j > 0) && (i != 1)) solver.mesh().connect_deformed(1, {serial_n[i][j-1], serial_n[i][j]}, {{1, 1}, {1, 0}});
        }
      }
    }
    int bc_sn = solver.mesh().add_boundary_condition(new hexed::Copy(), new Boundary_perturbation());
    solver.mesh().extrude();
    solver.mesh().cleanup();
    for (int i = 0; i < 3; ++i) solver.mesh().relax(.5); // so that we can see better
    solver.mesh().connect_rest(bc_sn);
    solver.calc_jacobian();
    hexed::Interpreter inter;
    solver.initialize(hexed::Spacetime_expr(hexed::Struct_expr("state0 = 0; state1 = 0; state2 = 1 + pos0; state3 = 1;"), inter));
    double area = solver.integral_field(hexed::Constant_func({1.}))[0];
    solver.mesh().disconnect_boundary(bc_sn);
    solver.mesh().extrude(true, .7);
    solver.mesh().cleanup();
    auto valid = solver.mesh().valid();
    REQUIRE(valid.n_redundant == 0);
    REQUIRE(valid.n_missing == 16);
    solver.mesh().connect_rest(bc_sn);
    solver.calc_jacobian();
    REQUIRE(solver.integral_field(Reciprocal_jacobian())[0] == Catch::Approx(40./4.)); // check number of elements and that jacobian is nonsingular
    REQUIRE(solver.integral_field(hexed::Constant_func({1.}))[0] == Catch::Approx(area)); // check that area has not changed since first extrusion
    // check that state variables have been interpolated correctly during second extrusion
    REQUIRE(solver.integral_field(hexed::Qpoint_expr(hexed::Struct_expr("errsq = (density - (1 + pos0))^2"), inter))[0] == Catch::Approx(0.).scale(1.));
  }
  SECTION("3D")
  {
    int serial_n [3][3][3];
    serial_n[1][1][1] = -1; // so that we know if we accidentally use this
    hexed::Solver solver {3, 2, 1.};
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
    solver.mesh().cleanup();
    auto valid = solver.mesh().valid();
    REQUIRE(valid.n_redundant == 0);
    REQUIRE(valid.n_missing == 60);
    int bc_sn = solver.mesh().add_boundary_condition(new hexed::Copy(), new hexed::Null_mbc());
    solver.mesh().connect_rest(bc_sn);
    solver.mesh().cleanup();
    solver.calc_jacobian();
    REQUIRE(solver.integral_field(Reciprocal_jacobian())[0] == Catch::Approx(86.)); // check number of elements
  }
}

TEST_CASE("normal continuity uncertainty")
{
  hexed::Solver sol({3, hexed::config::max_row_size, .2});
  int elem_sn = sol.mesh().add_element(0, true, {0, 0, 0});
  sol.mesh().extrude();
  int bc_sn = sol.mesh().add_boundary_condition(new hexed::Copy(), new hexed::Null_mbc());
  sol.mesh().connect_rest(bc_sn);
  sol.mesh().cleanup();
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  sol.set_uncert_surface_rep(bc_sn);
  for (auto handle : sol.mesh().elem_handles()) {
    double uncertainty = sol.sample(handle.ref_level, handle.is_deformed, handle.serial_n, hexed::Uncertainty())[0];
    if (handle.serial_n == elem_sn) REQUIRE(uncertainty == Catch::Approx(0).scale(1.));
    else CHECK(uncertainty == Catch::Approx(std::sqrt(2.)));
  }
}

TEST_CASE("artificial viscosity convergence")
{
  #if NDEBUG
  const int len0 = 100;
  const int len1 = 2;
  hexed::Solver sol(2, hexed::config::max_row_size, 1./len0);
  int sn [len0][len1];
  for (int i = 0; i < len0; ++i) {
    for (int j = 0; j < len1; ++j) {
      sn[i][j] = sol.mesh().add_element(0, 0, {i, j});
      if (i) sol.mesh().connect_cartesian(0, {sn[i - 1][j], sn[i][j]}, {0});
      if (j) sol.mesh().connect_cartesian(0, {sn[i][j - 1], sn[i][j]}, {1});
    }
  }
  int nonpen = sol.mesh().add_boundary_condition(new hexed::Nonpenetration, new hexed::Null_mbc);
  int pen [2];
  pen[0] = sol.mesh().add_boundary_condition(new hexed::Freestream(hexed::Mat<4>{1.1, 0., 1., 1.5}), new hexed::Null_mbc);
  pen[1] = sol.mesh().add_boundary_condition(new hexed::Freestream(hexed::Mat<4>{0.9, 0., 1., 1.5}), new hexed::Null_mbc);
  for (int positive = 0; positive < 2; ++positive) {
    for (int i = 0; i < len0; ++i) sol.mesh().connect_boundary(0, 0, sn[i][positive*(len1 - 1)], 1, positive, nonpen);
    for (int j = 0; j < len1; ++j) sol.mesh().connect_boundary(0, 0, sn[positive*(len0 - 1)][j], 0, positive, pen[positive]);
  }
  sol.mesh().cleanup();
  sol.mesh().valid().assert_valid();
  sol.calc_jacobian();
  double flow_width = 0.04;
  double adv_width = 0.01;
  sol.initialize(Tanh(flow_width));
  sol.nspace().assign("av_advect_iters", 6000);
  sol.nspace().assign("av_diff_iters", 3000);
  sol.nspace().assign("av_diff_ratio", 1e-6);
  sol.nspace().assign("freestream_speed", 1.);
  sol.nspace().assign("av_unscaled_max", 1000.);
  sol.nspace().assign("freestream" + std::to_string(2), 1.);
  sol.nspace().assign("freestream" + std::to_string(3), 1.5);
  sol.update_art_visc_smoothness(adv_width);
  REQUIRE(sol.iteration_status().adv_res < 1e-12);
  REQUIRE(sol.iteration_status().diff_res < 1e-12);
  double init_max = sol.bounds_field(hexed::Art_visc_coef())[0][1];
  // check that doubling the advection length multiplies the viscosity by ~2^max_row_size, with some wiggle room
  sol.update_art_visc_smoothness(adv_width*2);
  REQUIRE(sol.iteration_status().adv_res < 1e-12);
  REQUIRE(sol.iteration_status().diff_res < 1e-12);
  CHECK(std::log(sol.bounds_field(hexed::Art_visc_coef())[0][1]/init_max)/std::log(2) > hexed::config::max_row_size - 1.);
  // check that doubling the length scale of the flow divides the viscosity by ~2^(max_row_size - 1)
  sol.initialize(Tanh(flow_width/2));
  sol.update_art_visc_smoothness(adv_width);
  REQUIRE(sol.iteration_status().adv_res < 1e-12);
  REQUIRE(sol.iteration_status().diff_res < 1e-12);
  CHECK(std::log(sol.bounds_field(hexed::Art_visc_coef())[0][1]/init_max)/std::log(2) > hexed::config::max_row_size - 2.);
  #endif
}

TEST_CASE("uncertainty")
{
  hexed::Solver sol(2, 2, 1.);
  int sn = sol.mesh().add_element(0, true, {0, 0});
  sol.mesh().extrude();
  sol.mesh().extrude();
  int bc_sn = sol.mesh().add_boundary_condition(new hexed::Copy(), new hexed::Null_mbc());
  sol.mesh().connect_rest(bc_sn);
  sol.mesh().cleanup();
  sol.calc_jacobian();
  hexed::Position_func pos;
  sol.set_uncertainty(hexed::Elem_average(pos));
  // check that the uncertainty of the middle element is the x-coordinate of its centroid
  REQUIRE(sol.sample(0, true, sn, hexed::Uncertainty())[0] == Catch::Approx(.5));
  // uncertainty of the middle element will be set equal to the centroid of the (collapsed) element extruded 2 layers in the positive-x direction
  sol.synch_extruded_uncert();
  double correct = ((2*2*2 - 1.5*1.5*1.5)/3. - .5*(2.*2. - 1.5*1.5)/2.)/((1.5*1.5 - 1.)/2.);
  REQUIRE(sol.sample(0, true, sn, hexed::Uncertainty())[0] == Catch::Approx(correct));
}

TEST_CASE("cylinder tree mesh")
{
  static_assert(hexed::config::max_row_size >= 6);
  constexpr int row_size = 6;
  hexed::Solver solver (2, row_size, 1.);
  std::vector<hexed::Flow_bc*> bcs;
  hexed::Mat<2> origin{2., 2.};
  for (int i = 0; i < 4; ++i) bcs.push_back(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}));
  solver.mesh().add_tree(bcs, origin);
  for (int i = 0; i < 3; ++i) solver.mesh().update();
  solver.mesh().set_surface(new hexed::Hypersphere(origin, .5), new hexed::Nonpenetration, origin + Eigen::Vector2d{.8, .8});
  int n_initial = solver.mesh().n_elements();
  for (int i = 0; i < 3; ++i) solver.mesh().relax();
  solver.calc_jacobian();
  solver.initialize(hexed::Constant_func({0., 0., 1., 1e5}));
  REQUIRE_THAT(solver.integral_field(hexed::Constant_func({1.}))[0], Catch::Matchers::WithinRel(1 - M_PI*.25/4, 1e-6));
  for (int i = 0; i < 6; ++i) {
    // this criterion will refine all elements with a vertex that is within .1 of the midpoint of the arc
    auto criterion = [origin](hexed::Element& elem){
      bool ref = false;
      for (int i_vert = 0; i_vert < 4; ++i_vert) {
        double dist = 0;
        for (int i_dim = 0; i_dim < 2; ++i_dim) {
          dist += hexed::math::pow(elem.vertex(i_vert).pos[i_dim] - origin(i_dim) - .5/std::sqrt(2), 2);
        }
        double r = .1;
        ref = ref || dist < r*r;
      }
      ref = ref && elem.refinement_level() <= 6;
      return ref;
    };
    solver.mesh().update(criterion);
    for (int i = 0; i < 4; ++i) solver.mesh().relax();
    solver.mesh().valid().assert_valid();
  }
  solver.calc_jacobian();
  solver.initialize(hexed::Constant_func({0., 0., 1., 1e5}));
  REQUIRE_THAT(solver.integral_field(hexed::Constant_func({1.}))[0], Catch::Matchers::WithinRel(1 - M_PI*.25/4, 1e-6));
  for (int i = 0; i < 6; ++i) {
    solver.mesh().update(hexed::criteria::never, [](hexed::Element& elem){return elem.refinement_level() > 3;});
    for (int i = 0; i < 6; ++i) solver.mesh().relax();
    solver.mesh().valid().assert_valid();
  }
  REQUIRE(solver.mesh().n_elements() == n_initial); // this mesh should have been completely unrefined to where it started
  solver.calc_jacobian();
  solver.initialize(hexed::Constant_func({0., 0., 1., 1e5}));
  REQUIRE_THAT(solver.integral_field(hexed::Constant_func({1.}))[0], Catch::Matchers::WithinRel(1 - M_PI*.25/4, 1e-6));
}

TEST_CASE("file I/O")
{
  double tol = 1e-6;
  double correct_mass_integral;
  {
    hexed::Solver solver(2, hexed::config::max_row_size, .6);
    std::vector<hexed::Flow_bc*> bcs;
    for (int i = 0; i < 4; ++i) bcs.push_back(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}));
    solver.mesh().add_tree(bcs);
    solver.mesh().update();
    solver.mesh().set_surface(new hexed::Hypersphere(hexed::Mat<>::Zero(2), .2), new hexed::Nonpenetration, hexed::Mat<2>{.5, .5});
    solver.calc_jacobian();
    REQUIRE_THAT(solver.integral_field(hexed::Constant_func({1.}))[0], Catch::Matchers::WithinRel(.6*.6 - hexed::constants::pi*.2*.2/4, tol));
    solver.mesh().write("solver_io_test");
    solver.initialize(Nonuniform_mass());
    solver.write_state("solver_io_test");
    correct_mass_integral = solver.integral_field(hexed::Mass())[0];
  }
  SECTION("wrong params 0")
  {
    std::vector<hexed::Flow_bc*> bcs;
    for (int i = 0; i < 4; ++i) bcs.push_back(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}));
    hexed::Solver solver(1, hexed::config::max_row_size, .6);
    REQUIRE_THROWS(solver.read_mesh("solver_io_test", bcs, new hexed::Hypersphere(hexed::Mat<>::Zero(2), .2), new hexed::Nonpenetration));
  }
  SECTION("wrong params 1")
  {
    std::vector<hexed::Flow_bc*> bcs;
    for (int i = 0; i < 4; ++i) bcs.push_back(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}));
    hexed::Solver solver(2, hexed::config::max_row_size - 1, .6);
    REQUIRE_THROWS(solver.read_mesh("solver_io_test", bcs, new hexed::Hypersphere(hexed::Mat<>::Zero(2), .2), new hexed::Nonpenetration));
  }
  SECTION("right params")
  {
    std::vector<hexed::Flow_bc*> bcs;
    for (int i = 0; i < 4; ++i) bcs.push_back(new hexed::Freestream(Eigen::Vector4d{0., 0., 1., 1e5}));
    hexed::Solver solver(2, hexed::config::max_row_size, .6);
    solver.read_mesh("solver_io_test", bcs, new hexed::Hypersphere(hexed::Mat<>::Zero(2), .2), new hexed::Nonpenetration);
    REQUIRE_THAT(solver.integral_field(hexed::Constant_func({1.}))[0], Catch::Matchers::WithinRel(.6*.6 - hexed::constants::pi*.2*.2/4, tol));
    solver.read_state("solver_io_test");
    REQUIRE(solver.integral_field(hexed::Mass())[0] == Catch::Approx(correct_mass_integral));
  }
}
