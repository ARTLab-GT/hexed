#include <catch2/catch.hpp>

#include <cartdgConfig.hpp>
#include <Solution.hpp>
#include <get_gbc_cpg_euler.hpp>

class Supersonic_inlet : public cartdg::Ghost_boundary_condition
{
  public:
  int n_dim;

  Supersonic_inlet(cartdg::Grid& grid, int i_dim_arg, bool is_positive_face_arg)
  : cartdg::Ghost_boundary_condition(grid, i_dim_arg, is_positive_face_arg), n_dim(grid.n_dim)
  {}
    
  virtual void calc_ghost_state()
  {
    double mass = 2.;
    double speed = 700;
    double int_ener = 2e5;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim)
    {
      ghost_state().col(j_dim) = 0.;
    }
    ghost_state().col(i_dim) = mass*speed;
    if (is_positive_face) ghost_state().col(i_dim) *= -1;
    ghost_state().col(n_dim) = mass;
    ghost_state().col(n_dim + 1) = int_ener + 0.5*mass*speed*speed;
  }
};

TEST_CASE("Ghost boundary conditions")
{
  const int row_size = CARTDG_MAX_BASIS_ROW_SIZE;
  cartdg::Solution soln (5, 3, row_size, 1.);
  cartdg::Basis& basis = soln.basis;
  soln.kernel_settings.d_t_by_d_pos = 0.1;
  soln.add_block_grid(1, std::vector<int>{0, 0, 0}, std::vector<int>{3, 3, 3});
  cartdg::Grid& grid = soln.get_grid(0);
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_elem = 3*(k + 3*(j + 3*i));
        grid.pos[i_elem + 0] = i;
        grid.pos[i_elem + 1] = j;
        grid.pos[i_elem + 2] = k;
      }
    }
  }

  SECTION("Basic functionality")
  {
    const int n_face_qpoint = row_size*row_size;
    Supersonic_inlet bc (grid, 0, true);
    double some_number = 0.;
    REQUIRE(bc.jacobians.empty());
    REQUIRE(bc.elems.empty());
    bc.add_element(373, &some_number);
    bc.add_element(26);
    Supersonic_inlet new_bc (bc);
    for (cartdg::Ghost_boundary_condition* bc_ptr : {&bc, &new_bc})
    {
      cartdg::Ghost_boundary_condition& each_bc = *bc_ptr;
      REQUIRE(each_bc.i_dim == 0);
      REQUIRE(each_bc.n_var == 5);
      REQUIRE(each_bc.n_qpoint == n_face_qpoint);
      REQUIRE(each_bc.is_positive_face == true);
      REQUIRE(each_bc.domain_state().size() == 5*n_face_qpoint);
      REQUIRE(each_bc.ghost_state().size() == 5*n_face_qpoint);
      REQUIRE(each_bc.default_jacobian.size() == unsigned(3*3*grid.n_qpoint));
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        for (int i_dim : {0, 1, 2})
        {
          for (int j_dim : {0, 1, 2})
          {
            REQUIRE(each_bc.default_jacobian[(i_dim*3 + j_dim)*grid.n_qpoint + i_dim] == ((i_dim == j_dim) ? 1 : 0));
          }
        }
      }
      REQUIRE(each_bc.elems == std::vector<int> {373, 26});
      REQUIRE(each_bc.jacobians == std::vector<double*> {&some_number, bc.default_jacobian.data()});
    }
  }

  SECTION("Axis 0")
  {
    Supersonic_inlet bc0 (grid, 0, false);
    grid.ghost_bound_conds.push_back(&bc0);
    Supersonic_inlet bc1 (grid, 0, true);
    grid.ghost_bound_conds.push_back(&bc1);
    std::vector<double> other_jacobian (9*grid.n_qpoint, 0.);
    for (int i_dim = 0; i_dim < 3; ++i_dim)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        other_jacobian[i_dim*4*grid.n_qpoint + i_qpoint] = 3.;
      }
    }
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_elem; int i;
        i = 0;
        i_elem = k + 3*(j + 3*i);
        bc0.add_element(i_elem, other_jacobian.data());
        i = 2;
        i_elem = k + 3*(j + 3*i);
        bc1.add_element(i_elem);
      }
    }
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i = 0; i < row_size; ++i)
      {
        for (int j = 0; j < row_size; ++j)
        {
          for (int k = 0; k < row_size; ++k)
          {
            int i_qpoint = k + row_size*(j + row_size*i);
            #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
            #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
            READ(0) = 600.*(1 - 2.*i/(row_size - 1.)); READ(1) = 0.; READ(2) = 0.;
            READ(3) = 1.;
            READ(4) = 2e5;
            WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
            WRITE(3) = 1.;
            WRITE(4) = 2e5;
            #undef  READ
            #undef WRITE
          }
        }
      }
    }
    cartdg::get_gbc_cpg_euler(3, row_size)(grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                                       soln.basis, soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        if (std::abs(pos[i_qpoint]) < 1e-15)
        {
          REQUIRE(WRITE(3) == Approx(1. + 0.1*(700*2 - 600)/(3*basis.node_weights()[0])));
        }
        else if (std::abs(pos[i_qpoint] - 1.5) < 1e-15)
        {
          REQUIRE(WRITE(3) == Approx(1. + 0.1*(700*2 - 600)/basis.node_weights()[0]));
        }
        else
        {
          REQUIRE(WRITE(3) == 1.);
        }
        #undef WRITE
      }
    }
  }

  SECTION("Axis 1")
  {
    Supersonic_inlet bc1 (grid, 1, true);
    grid.ghost_bound_conds.push_back(&bc1);
    for (int i = 0; i < 3; ++i)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_elem; int j;
        j = 2;
        i_elem = k + 3*(j + 3*i);
        bc1.add_element(i_elem);
      }
    }
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        READ(0) = 0.; READ(1) = -600.;  READ(2) = 0.;
        READ(3) = 1.;
        READ(4) = 4e5;
        WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
        WRITE(3) = 1.;
        WRITE(4) = 4e5;
        #undef  READ
        #undef WRITE
      }
    }
    cartdg::get_gbc_cpg_euler(3, row_size)(grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                                       soln.basis, soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        double qpoint_pos = pos[i_qpoint + grid.n_qpoint];
        if (std::abs(qpoint_pos - 1.5) < 1e-15)
        {
          REQUIRE(WRITE(3) == Approx(1. + 0.1*(700*2 - 600)/basis.node_weights()[0]));
        }
        else
        {
          REQUIRE(WRITE(3) == 1.);
        }
        #undef WRITE
      }
    }
  }

  SECTION("Axis 2")
  {
    Supersonic_inlet bc0 (grid, 2, false);
    grid.ghost_bound_conds.push_back(&bc0);
    bc0.add_element(3);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        READ(0) = 0.; READ(1) = 0.; READ(2) = 600.;
        READ(3) = 1.;
        READ(4) = 4e5;
        WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
        WRITE(3) = 1.;
        WRITE(4) = 4e5;
        #undef  READ
        #undef WRITE
      }
    }
    cartdg::get_gbc_cpg_euler(3, row_size)(grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                                           soln.basis, soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        double qpoint_pos = pos[i_qpoint + 2*grid.n_qpoint];
        if ((i_elem == 3) && (std::abs(qpoint_pos) < 1e-15))
        {
          REQUIRE(WRITE(3) == Approx(1. + 0.1*(700*2 - 600)/basis.node_weights()[0]));
        }
        else
        {
          REQUIRE(WRITE(3) == 1.);
        }
        #undef WRITE
      }
    }
  }
}
