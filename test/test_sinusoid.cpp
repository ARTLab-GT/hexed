#include <cmath>

#include "catch.hpp"

#include "Grid.hpp"
#include "Equidistant.hpp"

TEST_CASE("Sinusoidal density wave")
{
  double length = 2*M_PI;
  double velocity = 9.;
  double pressure = 1e5;
  double mean_mass = 1.1;

  int n_divisions = 10;
  Equidistant basis (6);
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
        grid.state_r[i_qpoint + grid.n_dof*i_elem                  ] = mass;
        grid.state_r[i_qpoint + grid.n_dof*i_elem +   grid.n_qpoint] = momentum;
        grid.state_r[i_qpoint + grid.n_dof*i_elem + 2*grid.n_qpoint] = energy;
      }
    }
    //grid.print();
    grid.visualize("sinusoid_1d");
  }
}
