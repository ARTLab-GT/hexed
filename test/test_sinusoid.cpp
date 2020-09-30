#include <cmath>

#include <catch.hpp>

#include <Solution.hpp>
#include <Initializer.hpp>

class Sinusoidal_init : public Initializer
{
  public:
  const double length = 2*M_PI;
  const double pressure = 1e5;
  const double mean_mass = 1.1;
  const std::vector<double> velocity {0.5, 0.11, 0.7};
  int dim;

  Sinusoidal_init(int dim_arg) : dim(dim_arg) {}

  virtual std::vector<double> momentum(std::vector<double> position)
  {
    set_properties(position);
    return _momentum;
  }

  virtual std::vector<double> scalar_state(std::vector<double> position)
  {
    set_properties(position);
    std::vector<double> ss {mass, energy};
    return ss;
  }

  private:
  std::vector<double> _momentum;
  double mass;
  double energy;

  void set_properties(std::vector<double> pos)
  {
    double pos_func = 0.;
    for (int i = 0; i < dim; ++i)
    {
      pos_func += (i + 1)*pos[i];
    }
    mass = mean_mass + 0.1*sin(pos_func);
    energy = pressure/0.4;
    _momentum.clear();
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      _momentum.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
  }

};

TEST_CASE("Sinusoidal density wave")
{
  double length = 2.*M_PI;
  int rank = 6;
  double dt = 0.001;

  SECTION("1D")
  {
    Solution sol (3, 1, rank, length);
    sol.add_block_grid(2);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4};
    grid.auto_connect(periods);

    Sinusoidal_init init (1);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("sinusoid_1d");

    sol.update(dt);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.basis.rank; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*init.velocity[0]*cos(pos[i_qpoint]);
        double d_momentum_by_d_t = init.velocity[0]*d_mass_by_d_t;
        double d_energy_by_d_t = 0.5*init.velocity[0]*d_momentum_by_d_t;
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum_by_d_t
                == Approx(d_momentum_by_d_t).margin(0.01*init.mean_mass*init.velocity[0]));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*init.mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(init.pressure/0.4 +
                          0.5*init.mean_mass*init.velocity[0]*init.velocity[0])));
      }
    }
  }

  SECTION("2D")
  {
    Solution sol (4, 2, rank, length);
    sol.add_block_grid(3);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {8, 8};
    grid.auto_connect(periods);

    Sinusoidal_init init (2);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("sinusoid_2d");

    sol.update(dt);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*cos(pos[i_qpoint] + 2*pos[i_qpoint + grid.n_qpoint])
                               *(init.velocity[0] + 2*init.velocity[1]);
        double d_momentum0_by_d_t = init.velocity[0]*d_mass_by_d_t;
        double d_momentum1_by_d_t = init.velocity[1]*d_mass_by_d_t;
        double d_energy_by_d_t = 0.5*d_mass_by_d_t*(  init.velocity[0]*init.velocity[0]
                                                    + init.velocity[1]*init.velocity[1]);
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum0_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_momentum1_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum0_by_d_t
                == Approx(d_momentum0_by_d_t).margin(0.01*init.mean_mass*(  init.velocity[0]
                                                                          + init.velocity[1])));
        REQUIRE(num_d_momentum1_by_d_t
                == Approx(d_momentum1_by_d_t).margin(0.01*init.mean_mass*(  init.velocity[0]
                                                                          + init.velocity[1])));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*init.mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(init.pressure/0.4 +
                          0.5*init.mean_mass*(  init.velocity[0]*init.velocity[0]
                                              + init.velocity[1]*init.velocity[1]))));
      }
    }
  }

  SECTION("3D")
  {
    Solution sol (5, 3, rank, length);
    sol.add_block_grid(3);
    Grid& grid = sol.get_grid(0);
    std::vector<int> periods {8, 8, 8};
    grid.auto_connect(periods);

    Sinusoidal_init init (3);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    grid.visualize("sinusoid_3d");

    sol.update(dt);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        double d_mass_by_d_t = -0.1*cos(    pos[i_qpoint]
                                        + 2*pos[i_qpoint +   grid.n_qpoint]
                                        + 3*pos[i_qpoint + 2*grid.n_qpoint])
                               *(init.velocity[0] + 2*init.velocity[1] + 3*init.velocity[2]);
        double d_momentum0_by_d_t = init.velocity[0]*d_mass_by_d_t;
        double d_momentum1_by_d_t = init.velocity[1]*d_mass_by_d_t;
        double d_momentum2_by_d_t = init.velocity[2]*d_mass_by_d_t;
        double speed = sqrt(  init.velocity[0]*init.velocity[0]
                            + init.velocity[1]*init.velocity[1]
                            + init.velocity[2]*init.velocity[2]);
        double d_energy_by_d_t = 0.5*d_mass_by_d_t*speed*speed;
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum0_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_momentum1_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_momentum2_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum0_by_d_t
                == Approx(d_momentum0_by_d_t).margin(0.01*init.mean_mass*speed));
        REQUIRE(num_d_momentum1_by_d_t
                == Approx(d_momentum1_by_d_t).margin(0.01*init.mean_mass*speed));
        REQUIRE(num_d_momentum2_by_d_t
                == Approx(d_momentum2_by_d_t).margin(0.01*init.mean_mass*speed));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*init.mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(init.pressure/0.4 +
                          0.5*init.mean_mass*speed*speed)));
      }
    }
  }
}
