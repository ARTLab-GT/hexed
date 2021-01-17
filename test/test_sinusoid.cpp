#include <cmath>

#include <catch.hpp>

#include <Solution.hpp>

class Sinusoidal_init : public cartdg::Spacetime_func
{
  public:
  const double length = 2*M_PI;
  const double mean_pressure = 1e5;
  const double mean_mass = 1.1;
  const std::vector<double> velocity {0.5, 0.11, 0.7};
  int dim;

  Sinusoidal_init(int dim_arg) : dim(dim_arg) {}

  std::vector<double> operator()(std::vector<double> pos, double time)
  {
    pos_func = 0.;
    for (int i = 0; i < dim; ++i)
    {
      pos_func += (i + 1)*pos[i];
    }
    set_state();
    return state;
  }

  protected:
  double pos_func;
  std::vector<double> state;
  virtual void set_state() = 0;
};

class Velocity_init : public Sinusoidal_init
{
  public:
  Velocity_init(int dim_arg) : Sinusoidal_init(dim_arg) {}

  protected:
  virtual void set_state()
  {
    double mass = mean_mass + 0.1*sin(pos_func);
    double energy = mean_pressure/0.4;
    state.clear();
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      state.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
    state.push_back(mass);
    state.push_back(energy);
  }
};

class Pressure_init : public Sinusoidal_init
{
  public:
  Pressure_init(int dim_arg) : Sinusoidal_init(dim_arg) {}

  protected:
  virtual void set_state()
  {
    double mass = mean_mass;
    double pressure = mean_pressure*(1. + 0.01*sin(pos_func));
    double energy = pressure/0.4;
    state.clear();
    for (int i = 0; i < dim; ++i)
    {
      double momentum_i = velocity[i]*mass;
      state.push_back(momentum_i);
      energy += 0.5*momentum_i*velocity[i];
    }
    state.push_back(mass);
    state.push_back(energy);
  }
};

TEST_CASE("Sinusoidal density wave")
{
  double length = 2.*M_PI;
  int rank = 6;
  double cfl = 0.01;

  SECTION("1D")
  {
    cartdg::Solution sol (3, 1, rank, length);
    sol.add_block_grid(2);
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {4};
    grid.auto_connect(periods);

    Velocity_init init (1);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    sol.visualize("sinusoid_1d");

    double dt = sol.update(cfl);
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
                == Approx(d_energy_by_d_t).margin(0.01*(init.mean_pressure/0.4 +
                          0.5*init.mean_mass*init.velocity[0]*init.velocity[0])));
      }
    }

    Pressure_init pres_init (1);
    sol.initialize(pres_init);
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    sol.visualize("sinusoid_pressure_1d");

    dt = sol.update(cfl);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.basis.rank; ++i_qpoint)
      {
        double d_mass_by_d_t = 0.;
        double d_momentum_by_d_t = -0.01*init.mean_pressure*cos(pos[i_qpoint]);
        double d_energy_by_d_t = 0.5*init.velocity[0]*d_momentum_by_d_t;
        int i = i_qpoint + grid.n_dof*i_elem;
        double num_d_momentum_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_mass_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        double num_d_energy_by_d_t = (grid.state_r()[i] - initial[i])/dt;
        i += grid.n_qpoint;
        REQUIRE(num_d_momentum_by_d_t
                == Approx(d_momentum_by_d_t).margin(0.0001*init.mean_pressure));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*init.mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(init.mean_pressure/0.4 +
                          0.5*init.mean_mass*init.velocity[0]*init.velocity[0])));
      }
    }
  }

  SECTION("2D")
  {
    cartdg::Solution sol (4, 2, rank, length);
    sol.add_block_grid(3);
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {8, 8};
    grid.auto_connect(periods);

    Velocity_init init (2);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    sol.visualize("sinusoid_2d");

    double dt = sol.update(cfl);
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
                == Approx(d_energy_by_d_t).margin(0.01*(init.mean_pressure/0.4 +
                          0.5*init.mean_mass*(  init.velocity[0]*init.velocity[0]
                                              + init.velocity[1]*init.velocity[1]))));
      }
    }

    Pressure_init pres_init (2);
    sol.initialize(pres_init);
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    sol.visualize("sinusoid_2d");

    dt = sol.update(cfl);
    for (int i_elem = 0; i_elem < grid.n_elem; ++i_elem)
    {
      std::vector<double> pos = grid.get_pos(i_elem);
      for (int i_qpoint = 0; i_qpoint < grid.n_qpoint; ++i_qpoint)
      {
        double d_mass_by_d_t = 0.;
        double cos_pos = -0.01*cos(pos[i_qpoint] + 2*pos[i_qpoint + grid.n_qpoint])
                          *init.mean_pressure;
        double d_momentum0_by_d_t = cos_pos;
        double d_momentum1_by_d_t = cos_pos*2;
        double d_energy_by_d_t = 0.5*(  init.velocity[0]*d_momentum0_by_d_t
                                      + init.velocity[1]*d_momentum1_by_d_t);
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
                == Approx(d_momentum0_by_d_t).margin(0.0001*init.mean_pressure));
        REQUIRE(num_d_momentum1_by_d_t
                == Approx(d_momentum1_by_d_t).margin(0.0001*init.mean_pressure));
        REQUIRE(num_d_mass_by_d_t == Approx(d_mass_by_d_t).margin(0.01*init.mean_mass));
        REQUIRE(num_d_energy_by_d_t
                == Approx(d_energy_by_d_t).margin(0.01*(init.mean_pressure/0.4 +
                          0.5*init.mean_mass*(  init.velocity[0]*init.velocity[0]
                                              + init.velocity[1]*init.velocity[1]))));
      }
    }
  }

  SECTION("3D")
  {
    cartdg::Solution sol (5, 3, rank, length);
    sol.add_block_grid(3);
    cartdg::Grid& grid = sol.get_grid(0);
    std::vector<int> periods {8, 8, 8};
    grid.auto_connect(periods);

    Velocity_init init (3);
    sol.initialize(init);
    const int size = grid.n_elem*grid.n_dof;
    double initial[size];
    for (int i = 0; i < size; ++i)
    {
      initial[i] = grid.state_r()[i];
    }
    sol.visualize("sinusoid_3d");

    double dt = sol.update(cfl);
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
                == Approx(d_energy_by_d_t).margin(0.01*(init.mean_pressure/0.4 +
                          0.5*init.mean_mass*speed*speed)));
      }
    }
  }
}
