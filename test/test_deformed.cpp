#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Deformed_grid.hpp>
#include <Gauss_legendre.hpp>
#include <Tecplot_file.hpp>
#include <Ghost_boundary_condition.hpp>

class Copy_bc : public cartdg::Ghost_boundary_condition
{
  public:
  Copy_bc(cartdg::Storage_params params, int id, bool ip)
  : cartdg::Ghost_boundary_condition{params, id, ip}
  {}
  virtual void calc_ghost_state()
  {
    ghost_state() = domain_state();
  }
};

class Quadratic_def : public cartdg::Deformed_grid
{
  double deform (int i_node)
  {
    double node = basis.node(i_node);
    return node*(1. - node)*4.;
  }
  public:
  Quadratic_def(int n_var_arg, int n_dim_arg, int n_elem_arg, double mesh_size_arg, cartdg::Basis& basis_arg)
  : Deformed_grid(n_var_arg, n_dim_arg, n_elem_arg, mesh_size_arg, basis_arg) {}
  virtual std::vector<double> get_pos(int i_elem)
  {
    auto p {cartdg::Deformed_grid::get_pos(i_elem)};
    if (i_elem != 4) return p;
    double max_def [] {0.07, -0.13};
    for (int i_dim = 0; i_dim < 2; ++i_dim) {
      for (int i_row = 0; i_row < basis.row_size; ++i_row) {
        for (int j_row = 0; j_row < basis.row_size; ++j_row) {
          double def = max_def[i_dim]*mesh_size*deform(i_row)*deform(j_row);
          p[i_dim*n_qpoint + i_row*basis.row_size + j_row] += def;
        }
      }
    }
    return p;
  }
};

class Test
{
  protected:
  cartdg::Deformed_grid& grid;
  int row_size;
  int n_qpoint;
  int nfq;
  double veloc [2] {1.2, -1.3};
  double pres {101000};
  public:
  Test(cartdg::Deformed_grid& g) : grid{g}, row_size{grid.basis.row_size}, n_qpoint{row_size*row_size}, nfq{n_qpoint/row_size} {}

  void test()
  {
    for (int i = 0, i_elem = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        grid.add_element({i, j});
        ++i_elem;
      }
    }
    grid.element(4).time_step_scale()[0] = 0.77;
    connect();
    displace_vertices();
    for (int i = 0, i_elem = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        auto pos = grid.get_pos(i_elem);
        double* stage {grid.element(i_elem).stage(0)};
        for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint) {
          double mass = 1. + 0.03*pos[i_qpoint] + 0.02*pos[i_qpoint + n_qpoint];
          stage[0*n_qpoint + i_qpoint] = mass*veloc[0];
          stage[1*n_qpoint + i_qpoint] = mass*veloc[1];
          stage[2*n_qpoint + i_qpoint] = mass;
          stage[3*n_qpoint + i_qpoint] = pres/0.4 + 0.5*mass*(veloc[0]*veloc[0] + veloc[1]*veloc[1]);
        }
        ++i_elem;
      }
    }

    std::vector<Copy_bc> bcs;
    for (int i_dim : {0, 1}) {
      for (bool is_positive : {0, 1}) {
        bcs.emplace_back(grid.element(0).storage_params(), i_dim, is_positive);
      }
    }
    for (int i_dim : {0, 1}) {
      for (bool is_positive : {0, 1}) {
        for (int i_row = 0; i_row < 3; ++i_row) {
          int stride {i_dim ? 1 : 3};
          grid.add_element_gbc(i_row*3/stride + is_positive*2*stride, bcs[i_dim*2 + is_positive]);
        }
      }
    }
    grid.purge_vertices();
    grid.calc_jacobian();
    cartdg::Kernel_settings settings;
    double dt = 0.07;
    settings.d_t_by_d_pos = 0.07/0.2;

    grid.execute_write_face(settings);
    grid.execute_neighbor(settings);
    grid.execute_local(settings);
    for (int i_elem = 0; i_elem < 9; ++i_elem)
    {
      double* r {grid.element(i_elem).stage(settings.i_read)};
      double* w {grid.element(i_elem).stage(settings.i_write)};
      for (int i_dof = 0; i_dof < 4*n_qpoint; ++i_dof) {
        r[i_dof] = (w[i_dof] - r[i_dof])/dt;
      }
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint) {
        double tss = ((i_elem == 4) && (i_qpoint == 0)) ? 0.77 : 1.;
        REQUIRE(r[0*n_qpoint + i_qpoint] == Approx(tss*(-0.03*veloc[0]*veloc[0] - 0.02*veloc[1]*veloc[0])).scale(pres));
        REQUIRE(r[1*n_qpoint + i_qpoint] == Approx(tss*(-0.03*veloc[0]*veloc[1] - 0.02*veloc[1]*veloc[1])).scale(pres));
        REQUIRE(r[2*n_qpoint + i_qpoint] == Approx(tss*(-0.03*veloc[0] - 0.02*veloc[1])));
      }
    }
  }
  virtual void connect()
  {
    for (int i = 0, i_elem = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i > 0) grid.connect({i_elem - 3, i_elem}, {0, 0}, {1, 0});
        if (j > 0) grid.connect({i_elem - 1, i_elem}, {1, 1}, {1, 0});
        ++i_elem;
      }
    }
  }
  virtual void displace_vertices() {}
};

class Test_def_edge : public Test
{
  public:
  Test_def_edge(cartdg::Deformed_grid& g) : Test{g} {}
  virtual void displace_vertices()
  {
    grid.deformed_element(4).vertex(3).pos[0] += 0.05;
    grid.deformed_element(4).vertex(3).pos[1] -= 0.04;
  }
};

class Test_rotated : public Test_def_edge
{
  void connect(int i_elem, int i_dim)
  {
    std::array<int, 2> i_elems {i_elem - (3 - 2*i_dim), i_elem};
    std::array<int, 2> i_dims {i_dim, i_dim};
    std::array<bool, 2> is_p {1, 0};
    for (int i_side : {0, 1})
    {
      if (i_elems[i_side] == 4)
      {
        i_dims[i_side] = 1 - i_dims[i_side];
        if (i_dim) is_p[i_side] = !is_p[i_side];
      }
    }
    grid.connect(i_elems, i_dims, is_p);
  }
  public:
  Test_rotated(cartdg::Deformed_grid& g) : Test_def_edge{g} {}
  virtual void connect()
  {
    auto& elem {grid.deformed_element(4)};
    std::array<double, 3> pos {elem.vertex(0).pos};
    std::swap(elem.vertex(2).pos, pos);
    std::swap(elem.vertex(3).pos, pos);
    std::swap(elem.vertex(1).pos, pos);
    std::swap(elem.vertex(0).pos, pos);
    for (int i = 0, i_elem = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        if (i > 0) connect(i_elem, 0);
        if (j > 0) connect(i_elem, 1);
        ++i_elem;
      }
    }
  }
};

TEST_CASE("Deformed elements")
{
  const int row_size {cartdg::config::max_row_size};
  cartdg::Gauss_legendre basis {row_size};
  SECTION("plain")
  {
    cartdg::Deformed_grid grid {4, 2, 0, 0.2, basis};
    Test test {grid};
    test.test();
  }
  SECTION("deformed interior")
  {
    Quadratic_def grid {4, 2, 0, 0.2, basis};
    Test test {grid};
    test.test();
  }
  SECTION("deformed edges")
  {
    cartdg::Deformed_grid grid {4, 2, 0, 0.2, basis};
    Test_def_edge test {grid};
    test.test();
  }
  SECTION("rotated") // cyclical permutation of element 4's vertices
  {
    cartdg::Deformed_grid grid {4, 2, 0, 0.2, basis};
    Test_rotated test {grid};
    test.test();
  }
}
