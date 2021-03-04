#include <catch.hpp>

#include <Solution.hpp>
#include <kernels/neighbor/cpg_euler_gbc.hpp>

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
    double int_ener = 5e5;
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

TEST_CASE("Fitted boundary conditions")
{
  cartdg::Solution soln (5, 3, MAX_BASIS_RANK, 1.);
  cartdg::Basis& basis = soln.basis;
  soln.kernel_settings.d_t_by_d_pos = 1.0;
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

  SECTION("Axis 0")
  {
    Supersonic_inlet bc0 (grid, 0, false);
    grid.ghost_bound_conds.push_back(&bc0);
    Supersonic_inlet bc1 (grid, 0, true);
    grid.ghost_bound_conds.push_back(&bc1);
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        int i_elem; int i;
        i = 0;
        i_elem = k + 3*(j + 3*i);
        bc0.elems.push_back(i_elem);
        i = 2;
        i_elem = k + 3*(j + 3*i);
        bc1.elems.push_back(i_elem);
      }
    }
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        READ(0) = 0.; READ(1) = 0.; READ(2) = 0.;
        READ(3) = 1.;
        READ(4) = 4e5;
        WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
        WRITE(3) = 1.;
        WRITE(4) = 4e5;
        #undef  READ
        #undef WRITE
      }
    }
    auto weights = soln.basis.node_weights();
    cartdg::cpg_euler_gbc<5, MAX_BASIS_RANK*MAX_BASIS_RANK*MAX_BASIS_RANK, MAX_BASIS_RANK>
                         (grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                         weights[0], soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        if (std::min(std::abs(pos[i_qpoint]), std::abs(pos[i_qpoint] - 1.5)) < 1e-15)
        {
          REQUIRE(WRITE(3) == 1. + 2.*700/basis.node_weights()[0]);
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
        bc1.elems.push_back(i_elem);
      }
    }
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        READ(0) = 0.; READ(1) = 0.; READ(2) = 0.;
        READ(3) = 1.;
        READ(4) = 4e5;
        WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
        WRITE(3) = 1.;
        WRITE(4) = 4e5;
        #undef  READ
        #undef WRITE
      }
    }
    auto weights = soln.basis.node_weights();
    cartdg::cpg_euler_gbc<5, MAX_BASIS_RANK*MAX_BASIS_RANK*MAX_BASIS_RANK, MAX_BASIS_RANK>
                         (grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                         weights[0], soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        double qpoint_pos = pos[i_qpoint + grid.n_qpoint];
        if (std::abs(qpoint_pos - 1.5) < 1e-15)
        {
          REQUIRE(WRITE(3) == 1. + 2.*700/basis.node_weights()[0]);
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
    bc0.elems.push_back(3);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define  READ(i) grid.state_r()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        READ(0) = 0.; READ(1) = 0.; READ(2) = 0.;
        READ(3) = 1.;
        READ(4) = 4e5;
        WRITE(0) = 0.; WRITE(1) = 0.; WRITE(2) = 0.;
        WRITE(3) = 1.;
        WRITE(4) = 4e5;
        #undef  READ
        #undef WRITE
      }
    }
    auto weights = soln.basis.node_weights();
    cartdg::cpg_euler_gbc<5, MAX_BASIS_RANK*MAX_BASIS_RANK*MAX_BASIS_RANK, MAX_BASIS_RANK>
                         (grid.ghost_bound_conds, grid.state_r(), grid.state_w(),
                          weights[0], soln.kernel_settings);
    for (int i_elem = 0; i_elem < 27; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        #define WRITE(i) grid.state_w()[i_elem*grid.n_dof + i_qpoint + (i)*grid.n_qpoint]
        double qpoint_pos = pos[i_qpoint + 2*grid.n_qpoint];
        if ((i_elem == 3) && (std::abs(qpoint_pos) < 1e-15))
        {
          REQUIRE(WRITE(3) == 1. + 2.*700/basis.node_weights()[0]);
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
