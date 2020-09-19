#include <cmath>

#include "catch.hpp"

#include "Grid.hpp"
#include "Equidistant.hpp"
#include "kernels/local/cpg_euler_matrix.hpp"

TEST_CASE("Sinusoidal density wave")
{
  double length = 2*M_PI;
  double velocity = 9.;
  double pressure = 1e5;
  double mean_mass = 1.1;

  int n_divisions = 5;
  const int rank = 6;
  Equidistant basis (rank);
  SECTION("1D")
  {
    Grid grid (3, 1, n_divisions, length/n_divisions, basis);
    for (int i_elem = 0; i_elem < n_divisions; ++i_elem)
    {
      grid.pos[i_elem] = i_elem;
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
      {
        double mass = mean_mass + 0.1*sin(pos[i_qpoint]);
        double momentum = velocity*mass;
        double energy = pressure/0.4 + 0.5*momentum*velocity;
        grid.state_r[i_qpoint + grid.n_dof*i_elem                  ] = momentum;
        grid.state_r[i_qpoint + grid.n_dof*i_elem +   grid.n_qpoint] = mass;
        grid.state_r[i_qpoint + grid.n_dof*i_elem + 2*grid.n_qpoint] = energy;
      }
    }
    //grid.print();
    grid.visualize("sinusoid_1d_0");
    cpg_euler_matrix<3, rank, rank>(grid.state_r, grid.state_w, n_divisions,
                                    grid.basis.diff_mat(), n_divisions/length);
    for (int i_elem = 0; i_elem < n_divisions; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*velocity*cos(pos[i_qpoint]);
        double d_momentum_by_d_t = velocity*d_mass_by_d_t;
        double d_energy_by_d_t = 0.5*velocity*d_momentum_by_d_t;
        double num_d_mass_by_d_t =   grid.state_w[i_qpoint + grid.n_dof*i_elem
                                                  + grid.n_qpoint]
                                   - grid.state_r[i_qpoint + grid.n_dof*i_elem
                                                  + grid.n_qpoint];
        grid.state_r[i_qpoint + grid.n_dof*i_elem + grid.n_qpoint] = num_d_mass_by_d_t;
        double num_d_momentum_by_d_t =   grid.state_w[i_qpoint + grid.n_dof*i_elem]
                                       - grid.state_r[i_qpoint + grid.n_dof*i_elem];
        double num_d_energy_by_d_t =   grid.state_w[i_qpoint + grid.n_dof*i_elem
                                                   + 2*grid.n_qpoint]
                                    - grid.state_r[i_qpoint + grid.n_dof*i_elem
                                                   + 2*grid.n_qpoint];
        CHECK(num_d_momentum_by_d_t
                == Approx(d_momentum_by_d_t).margin(0.01*mean_mass*velocity));
        CHECK(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*mean_mass));
        CHECK(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(pressure/0.4 +
                          0.5*mean_mass*velocity*velocity)));
      }
    }
  }
}
