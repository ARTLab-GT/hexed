#include <cmath>
#include <iostream>

#include <catch.hpp>

#include <Grid.hpp>
#include <Equidistant.hpp>
#include <kernels/local/cpg_euler_matrix.hpp>

TEST_CASE("Sinusoidal density wave")
{
  double length = 2*M_PI;
  double pressure = 1e5;
  double mean_mass = 1.1;

  int n_divisions = 5;
  const int rank = 6;
  Equidistant basis (rank);
  SECTION("1D")
  {
    double velocity = 9.;
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
        grid.state_r()[i_qpoint + grid.n_dof*i_elem                  ] = momentum;
        grid.state_r()[i_qpoint + grid.n_dof*i_elem +   grid.n_qpoint] = mass;
        grid.state_r()[i_qpoint + grid.n_dof*i_elem + 2*grid.n_qpoint] = energy;
      }
    }
    //grid.print();
    grid.visualize("sinusoid_1d");
    cpg_euler_matrix<3, rank, rank>(grid.state_r(), grid.state_w(), n_divisions,
                                    grid.basis.diff_mat(), n_divisions/length);
    for (int i_elem = 0; i_elem < n_divisions; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < basis.rank; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*velocity*cos(pos[i_qpoint]);
        double d_momentum_by_d_t = velocity*d_mass_by_d_t;
        double d_energy_by_d_t = 0.5*velocity*d_momentum_by_d_t;
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum_by_d_t
                == Approx(d_momentum_by_d_t).margin(0.01*mean_mass*velocity));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(pressure/0.4 +
                          0.5*mean_mass*velocity*velocity)));
      }
    }
  }

  SECTION("2D")
  {
    double velocity0 = 0.3;
    double velocity1 = 0.7;
    const int n_qpoint = rank*rank;
    int n_elem = n_divisions*n_divisions;
    Grid grid (4, 2, n_elem, length/n_divisions, basis);
    for (int i_row = 0; i_row < n_divisions; ++i_row)
    {
      for (int i_col = 0; i_col < n_divisions; ++i_col)
      {
        int i_elem = i_row*n_divisions + i_col;
        grid.pos[2*i_elem] = i_col;
        grid.pos[2*i_elem + 1] = i_row;
        std::vector<double> pos = grid.get_pos(i_elem);
        for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
        {
          double mass = mean_mass + 0.1*sin(2*pos[i_qpoint] + pos[i_qpoint + n_qpoint]);
          double momentum0 = velocity0*mass;
          double momentum1 = velocity1*mass;
          double energy = pressure/0.4 + 0.5*(momentum0*velocity0 + momentum1*velocity1);
          grid.state_r()[i_qpoint + grid.n_dof*i_elem                  ] = momentum0;
          grid.state_r()[i_qpoint + grid.n_dof*i_elem +   grid.n_qpoint] = momentum1;
          grid.state_r()[i_qpoint + grid.n_dof*i_elem + 2*grid.n_qpoint] = mass;
          grid.state_r()[i_qpoint + grid.n_dof*i_elem + 3*grid.n_qpoint] = energy;
        }
      }
    }
    //grid.print();
    grid.visualize("sinusoid_2d");
    cpg_euler_matrix<4, n_qpoint, rank>(grid.state_r(), grid.state_w(), n_elem,
                                        grid.basis.diff_mat(), n_divisions/length);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*cos(2*pos[i_qpoint] + pos[i_qpoint + n_qpoint])
                               *(2*velocity0 + velocity1);
        double d_momentum0_by_d_t = velocity0*d_mass_by_d_t;
        double d_momentum1_by_d_t = velocity1*d_mass_by_d_t;
        double d_energy_by_d_t = 0.5*d_mass_by_d_t*(velocity0*velocity0 + velocity1*velocity1);
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum0_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_momentum1_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum0_by_d_t
                == Approx(d_momentum0_by_d_t).margin(0.01*mean_mass*(velocity0 + velocity1)));
        REQUIRE(num_d_momentum1_by_d_t
                == Approx(d_momentum1_by_d_t).margin(0.01*mean_mass*(velocity0 + velocity1)));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(pressure/0.4 +
                          0.5*mean_mass*(velocity0*velocity0 + velocity1*velocity1))));
      }
    }
  }

  SECTION("3D")
  {
    double velocity0 = 0.5;
    double velocity1 = 0.11;
    double velocity2 = 0.7;
    const int n_qpoint = rank*rank*rank;
    int n_elem = n_divisions*n_divisions*n_divisions;
    Grid grid (5, 3, n_elem, length/n_divisions, basis);
    for (int i_row = 0; i_row < n_divisions; ++i_row)
    {
      for (int j_row = 0; j_row < n_divisions; ++j_row)
      {
        for (int k_row = 0; k_row < n_divisions; ++k_row)
        {
          int i_elem = (i_row*n_divisions + j_row)*n_divisions + k_row;
          grid.pos[3*i_elem    ] = i_row;
          grid.pos[3*i_elem + 1] = j_row;
          grid.pos[3*i_elem + 2] = k_row;
          std::vector<double> pos = grid.get_pos(i_elem);
          for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
          {
            double mass = mean_mass + 0.1*sin(    pos[i_qpoint]
                                              + 2*pos[i_qpoint +   n_qpoint]
                                              + 3*pos[i_qpoint + 2*n_qpoint]);
            double momentum0 = velocity0*mass;
            double momentum1 = velocity1*mass;
            double momentum2 = velocity2*mass;
            double energy = pressure/0.4 + 0.5*(+ momentum0*velocity0
                                                + momentum1*velocity1
                                                + momentum2*velocity2);
            grid.state_r()[i_qpoint + grid.n_dof*i_elem                  ] = momentum0;
            grid.state_r()[i_qpoint + grid.n_dof*i_elem +   grid.n_qpoint] = momentum1;
            grid.state_r()[i_qpoint + grid.n_dof*i_elem + 2*grid.n_qpoint] = momentum2;
            grid.state_r()[i_qpoint + grid.n_dof*i_elem + 3*grid.n_qpoint] = mass;
            grid.state_r()[i_qpoint + grid.n_dof*i_elem + 4*grid.n_qpoint] = energy;
          }
        }
      }
    }
    //grid.print();
    grid.visualize("sinusoid_3d");
    cpg_euler_matrix<5, n_qpoint, rank>(grid.state_r(), grid.state_w(), n_elem,
                                        grid.basis.diff_mat(), n_divisions/length);
    for (int i_elem = 0; i_elem < n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < n_qpoint; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*cos(    pos[i_qpoint]
                                        + 2*pos[i_qpoint +   n_qpoint]
                                        + 3*pos[i_qpoint + 2*n_qpoint])
                               *(velocity0 + 2*velocity1 + 3*velocity2);
        double d_momentum0_by_d_t = velocity0*d_mass_by_d_t;
        double d_momentum1_by_d_t = velocity1*d_mass_by_d_t;
        double d_momentum2_by_d_t = velocity2*d_mass_by_d_t;
        double speed = sqrt(velocity0*velocity0 + velocity1*velocity1 + velocity2*velocity2);
        double d_energy_by_d_t = 0.5*d_mass_by_d_t*speed*speed;
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum0_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_momentum1_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_momentum2_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = grid.state_w()[i] - grid.state_r()[i];
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum0_by_d_t
                == Approx(d_momentum0_by_d_t).margin(0.01*mean_mass*speed));
        REQUIRE(num_d_momentum1_by_d_t
                == Approx(d_momentum1_by_d_t).margin(0.01*mean_mass*speed));
        REQUIRE(num_d_momentum2_by_d_t
                == Approx(d_momentum2_by_d_t).margin(0.01*mean_mass*speed));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(pressure/0.4 +
                          0.5*mean_mass*speed*speed)));
      }
    }
  }
}
