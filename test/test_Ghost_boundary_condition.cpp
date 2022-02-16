#include <catch2/catch.hpp>

#include <Ghost_boundary_condition.hpp>
#include <cartdgConfig.hpp>
#include <Solution.hpp>
#include <get_gbc_convective.hpp>
#include <get_gbc_deformed_convective.hpp>

class Supersonic_inlet : public cartdg::Ghost_boundary_condition
{
  public:
  int n_dim;

  Supersonic_inlet(cartdg::Storage_params params, int i_dim_arg, bool is_positive_face_arg)
  : cartdg::Ghost_boundary_condition(params, i_dim_arg, is_positive_face_arg), n_dim(params.n_dim)
  {}

  virtual void calc_ghost_state()
  {
    double mass = 2.;
    double speed = 700;
    double int_ener = 2e5;
    for (int j_dim = 0; j_dim < n_dim; ++j_dim) {
      ghost_state().col(j_dim) = 0.;
    }
    ghost_state().col(i_dim()) = mass*speed;
    if (is_positive_face()) ghost_state().col(i_dim()) *= -1;
    ghost_state().col(n_dim) = mass;
    ghost_state().col(n_dim + 1) = int_ener + 0.5*mass*speed*speed;
  }
};

TEST_CASE("Ghost_boundary_condition class")
{
  const int row_size = cartdg::config::max_row_size;
  const int n_face_qpoint = row_size*row_size;
  cartdg::Gauss_legendre basis {row_size};
  cartdg::Regular_grid grid {5, 3, 0, 1., basis};
  cartdg::Storage_params params {3, 5, 3, row_size};
  Supersonic_inlet bc (params, 0, true);
  Supersonic_inlet new_bc (bc);
  for (cartdg::Ghost_boundary_condition* bc_ptr : {&bc, &new_bc}) {
    cartdg::Ghost_boundary_condition& each_bc = *bc_ptr;
    REQUIRE(each_bc.domain_state().size() == 5*n_face_qpoint);
    REQUIRE(each_bc.ghost_state().size() == 5*n_face_qpoint);
  }
}

TEST_CASE("gbc kernels")
{
  const int row_size = cartdg::config::max_row_size;
  cartdg::Solution soln (5, 3, row_size, 1.);
  soln.kernel_settings.d_t_by_d_pos = 0.1;
  soln.add_block_grid(1, std::vector<int>{0, 0, 0}, std::vector<int>{3, 3, 3});
  cartdg::Regular_grid& grid = soln.reg_grids[0];
  soln.add_deformed_grid(1.);
  cartdg::Deformed_grid& def_grid = soln.def_grids[0];
  def_grid.add_element({0, 0, 0});
  def_grid.deformed_element(0).vertex(1).pos[2] /= 2.; // make the face half as large
  def_grid.deformed_element(0).vertex(3).pos[2] /= 2.;

  cartdg::Storage_params params {3, 5, 3, row_size};
  Supersonic_inlet gbc0 {params, 0, false};
  Supersonic_inlet gbc1 {params, 0, true};
  Supersonic_inlet gbc2 {params, 2, false};
  grid.add_element_gbc(0, gbc0);
  grid.add_element_gbc(1, gbc0);
  grid.add_element_gbc(26, gbc1);
  grid.add_element_gbc(0, gbc2);
  def_grid.add_element_gbc(0, gbc0);
  def_grid.calc_jacobian();

  int nfqpoint = row_size*row_size;
  for (cartdg::Grid* g : soln.all_grids()) {
    for (int i_elem = 0; i_elem < g->n_elem; ++i_elem) {
      double* face = g->element(i_elem).face();
      for (int i_dim : {0, 1, 2}) {
        for (int positive : {0, 1}) {
          for (int i_qpoint = 0; i_qpoint < nfqpoint; ++i_qpoint) {
            double* qpoint = face + (2*i_dim + positive)*5*nfqpoint + i_qpoint;
            qpoint[0*nfqpoint] = 700.;
            qpoint[1*nfqpoint] = 0.;
            qpoint[2*nfqpoint] = 700.;
            qpoint[3*nfqpoint] = 1.;
            qpoint[4*nfqpoint] = 1e5 + 0.5*1.*2*700*700;
          }
        }
      }
    }
  }

  grid.execute_neighbor(soln.kernel_settings);
  REQUIRE(grid.element(0).face()[3*nfqpoint] == Approx(2.*700.));
  REQUIRE(grid.element(0).face()[3*nfqpoint + nfqpoint - 1] == Approx(2.*700.));
  REQUIRE(grid.element(1).face()[3*nfqpoint] == Approx(2.*700.));
  REQUIRE(grid.element(26).face()[(1*5 + 3)*nfqpoint] < 0.);
  REQUIRE(grid.element(0).face()[(2*2*5 + 3)*nfqpoint] == Approx(2.*700.));
  REQUIRE(def_grid.element(0).face()[3*nfqpoint] == Approx(2.*700./2.)); // flux should be proportional to face size
}
