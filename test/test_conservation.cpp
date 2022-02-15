#include <sys/stat.h>
#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Solution.hpp>

class Initializer : public cartdg::Spacetime_func
{
  public:
  int dim;
  const std::vector<double> velocity {4., -.3, 1.2};
  Initializer(int dim_arg) : dim(dim_arg) {}

  std::vector<double> operator()(std::vector<double> position, double time)
  {
    double mass = 1.;
    double pressure = 1.e5;
    double energy = pressure/0.4;
    std::vector<double> state;
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      state.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
    state.push_back(mass);
    state.push_back(energy);
    return state;
  }
};

TEST_CASE("Conservation of state variables")
{
  double length = 1.;
  int row_size = 3;

  SECTION("1D")
  {
    cartdg::Solution sol (3, 1, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4};
    grid.auto_connect(periods);

    Initializer init (1);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D cartesian")
  {
    cartdg::Solution sol (4, 2, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4, 4};
    grid.auto_connect(periods);

    Initializer init (2);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D deformed")
  {
    const int row_size = cartdg::config::max_row_size;
    cartdg::Solution sol (4, 2, row_size, 1.);
    sol.add_empty_grid(1);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    sol.add_deformed_grid(1);
    cartdg::Deformed_grid& def_grid = sol.def_grids[0];

    // I know this setup is complicated... sorry :(
    /*
     * Element diagram.
     * [i] = grid.element(i)
     * /i/ = def_grid.deformed_element(i)
     *
     * x1  ^
     * 2.5 _| [ 5][16][17][18][19][11]
     * 2.0 _| [ 4]/ 3// 7//11//15/[10]
     * 1.5 _| [ 3]/ 2// 6//10//14/[ 9]
     * 1.0 _| [ 2]/ 1// 5// 9//13/[ 8]
     * 0.5 _| [ 1]/ 0// 4// 8//12/[ 7]
     * 0.0 _| [ 0][12][13][14][15][ 6]
     *      + ------------------------>
     *        |   |   |   |   |   |
     *        0.0 0.5 1.0 1.5 2.0 2.5 x0
     */
    for (int i = 0; i < 6; ++i) grid.add_element({0, i});
    for (int i = 0; i < 6; ++i) grid.add_element({5, i});
    for (int i = 0; i < 4; ++i) grid.add_element({i + 1, 0});
    for (int i = 0; i < 4; ++i) grid.add_element({i + 1, 5});
    grid.auto_connect({6, 6});

    for (int i = 1; i < 5; ++i) {
      for (int j = 1; j < 5; ++j) {
        def_grid.add_element({i, j});
      }
    }
    double center [] {1.5 + .1, 1.5 + .1};
    {
      // move one vertex to create deformation
      cartdg::Deformed_element& elem {def_grid.deformed_element(5)};
      elem.vertex(3).pos = {center[0], center[1], 0.};
      elem.node_adjustments()[row_size + 1] = 0.1;
    }
    // deform corresponding vertices to match
    {
      cartdg::Deformed_element& elem {def_grid.deformed_element(9)};
      elem.vertex(1).pos = {center[0], center[1], 0.};
    }
    {
      cartdg::Deformed_element& elem {def_grid.deformed_element(10)};
      elem.vertex(0).pos = {center[0], center[1], 0.};
    }
    {
      // rotate element to facilitate testing connections along different axes
      cartdg::Deformed_element& elem {def_grid.deformed_element(6)};
      elem.vertex(0).pos = {1.0, 2.0, 0.};
      elem.vertex(1).pos = {1.5, 2.0, 0.};
      elem.vertex(2).pos = {1.0, 1.5, 0.};
      elem.vertex(3).pos = {center[0], center[1], 0.};
    }
    {
      // perturb another vertex (don't bother to move corresponding vertices)
      cartdg::Deformed_element& elem {def_grid.deformed_element(15)};
      elem.vertex(0).pos[0] += 0.2;
      elem.vertex(0).pos[1] += 0.2;
    }
    // connect all the deformed elements except 6 (the rotated one)
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
        int inds [2][2] {{4*i+j, 4*(i+1)+j}, {4*j+i, 4*j+i+1}};
        for (int i_dim : {0, 1}) {
          if ((inds[i_dim][0] != 6) && (inds[i_dim][1] != 6)) {
            def_grid.connect({inds[i_dim][0], inds[i_dim][1]}, {i_dim, i_dim}, {1, 0});
          }
        }
      }
    }
    // connect element 6
    def_grid.connect({6, 7}, {0, 1}, {0, 0});
    def_grid.connect({6, 5}, {0, 1}, {1, 1});
    def_grid.connect({6,  2}, {1, 0}, {0, 1});
    def_grid.connect({6, 10}, {1, 0}, {1, 0});
    // connect deformed to non-deformed
    for (int i = 0; i < 4; ++i) {
      grid.add_connection(&grid.element(i + 1), &def_grid.element(i), 0);
      grid.add_connection(&def_grid.element(i + 12), &grid.element(i + 7), 0);
      grid.add_connection(&grid.element(i + 12), &def_grid.element(4*i), 1);
      grid.add_connection(&def_grid.element(4*i + 3), &grid.element(i + 16), 1);
    }
    def_grid.calc_jacobian();

    Initializer init (2);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state) {
      for (int i_elem : {5, 6, 15}) def_grid.element(i_elem).stage(0)[i_state] = i_state + 1;
    }

    auto before {sol.integral()};
    sol.visualize_field("deformed_conservation");
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("2D hanging node")
  {
    cartdg::Solution sol {4, 2, row_size, length};
    sol.add_empty_grid(1);
    sol.add_empty_grid(2);
    cartdg::Regular_grid& grid1 {sol.reg_grids[0]};
    cartdg::Regular_grid& grid2 {sol.reg_grids[1]};
    grid1.add_element({1, 1});
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        if ((i < 2) || (j < 2)) grid2.add_element({i, j});
      }
    }
    grid2.auto_connect({4, 4});
    grid2.connect_refined(&grid1.element(0), {&grid2.element(2), &grid2.element( 3)}, 0, 1);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(6), &grid2.element( 7)}, 0, 0);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(8), &grid2.element(10)}, 1, 1);
    grid2.connect_refined(&grid1.element(0), {&grid2.element(9), &grid2.element(11)}, 1, 0);
    Initializer init {2};
    sol.initialize(init);
    for (cartdg::Regular_grid& grid : sol.reg_grids) {
      for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem) {
        for (int i_dof = 0; i_dof < grid.n_dof; ++i_dof) {
          grid.element(i_elem).stage(0)[i_dof] += 0.01*(std::rand()%10);
        }
      }
    }
    auto before {sol.integral()};
    sol.visualize_field("hanging_node");
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < sol.n_var; ++i_var) {
      REQUIRE((before[i_var] - after[i_var])/dt == Approx(0).scale(std::abs(before[i_var]) + 1.));
    }
  }

  SECTION("3D")
  {
    cartdg::Solution sol (5, 3, row_size, length);
    sol.add_block_grid(2);
    cartdg::Regular_grid& grid = sol.reg_grids[0];
    std::vector<int> periods {4, 4, 4};
    grid.auto_connect(periods);

    Initializer init (3);
    sol.initialize(init);
    for (int i_state = 0; i_state < grid.n_qpoint; ++i_state)
    {
      grid.element(0).stage(0)[i_state] = i_state + 1;
    }
    auto before {sol.integral()};
    double dt = sol.update();
    auto after {sol.integral()};
    for (int i_var = 0; i_var < grid.n_var; ++i_var) {
      REQUIRE((after[i_var] - before[i_var])/dt == Approx{0.}.scale(std::abs(before[i_var]) + 1.));
    }
  }
}
