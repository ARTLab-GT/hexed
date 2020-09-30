#include <iostream>
#include <chrono>

#include <Initializer.hpp>
#include <Solution.hpp>

class Vortex_init : public Initializer
{
  public:
  // Physical parameters
  double sp_heat_rat = 1.4;
  double sp_gas_const = 287;

  // Freestream parameters
  double freestream_mass = 1.;
  double freestream_velocity = 100.;
  double freestream_pressure = 1.e5;

  // Vortex parameters
  double vortex_strength = 30;
  double critical_radius = 0.05;
  double decay_rate = 0.2;

  // Geometric parameters --
  // if you mess with these you must also modify the Solution setup in main()
  const double center0 = 0.5;
  const double center1 = 0.5;
  const int dim = 2;

  // Can use this for estimating cfl number
  double max_char_speed = 0.;

  Vortex_init() { _momentum.resize(2); }

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

  double sqr(double x)
  {
    return x*x;
  }

  void set_properties(std::vector<double> pos)
  {
    // Calculate vortex perturbation
    double radius = sqrt(sqr(pos[0] - center0) + sqr(pos[1] - center1));
    double nondim_radius = radius/critical_radius;
    double veloc_mag = vortex_strength*nondim_radius*exp(decay_rate*(1. - sqr(nondim_radius)));
    double veloc0 =  veloc_mag*(pos[1] - center1 + 1e-6)/(radius + 1e-6);
    double veloc1 = -veloc_mag*(pos[0] - center0 + 1e-6)/(radius + 1e-6);
    double temperature = -(sp_heat_rat - 1)/(4*decay_rate*sp_heat_rat*sp_gas_const)
                         *sqr(vortex_strength*exp(decay_rate*(1. - sqr(nondim_radius))));

    // Combine with freestream
    veloc0 += freestream_velocity;
    double freestream_temperature = freestream_pressure/(sp_gas_const*freestream_mass);
    temperature += freestream_temperature;

    // Calculate char speed
    double char_speed = veloc_mag + sqrt(sp_heat_rat*sp_gas_const*temperature);
    max_char_speed = std::max(max_char_speed, char_speed);

    // Convert to conserved variables
    mass = freestream_mass*pow(temperature/freestream_temperature, 1./(sp_heat_rat - 1.));
    _momentum[0] = mass*veloc0; _momentum[1] = mass*veloc1;
    energy = sp_gas_const/(sp_heat_rat - 1.)*temperature + 0.5*mass*(sqr(veloc0) + sqr(veloc1));
  }

};

int main()
{
  // Resolution parameters
  const int rank = 6;
  const int ref_level = 3;

  // Solution setup
  Solution solution (4, 2, rank, 1.);
  solution.add_block_grid(ref_level);
  Grid& grid = solution.get_grid(0);
  int n_div = 1; for (int i = 0; i < ref_level; ++i) n_div *= 2;
  std::vector<int> periods {n_div, n_div};
  grid.auto_connect(periods);
  Vortex_init v_init;
  solution.initialize(v_init);

  // Let's go!
  grid.visualize("initial");
  auto start = std::chrono::high_resolution_clock::now();
  double time = 0;
  for (int i = 0; i < 20; ++i)
  {
    time += 0.001;
    while (grid.time < time)
    {
      solution.update();
    }
    char buffer [100];
    snprintf(buffer, 100, "t%e.txt", time);
    grid.visualize(std::string(buffer));
  }
  time = 0.2;
  while (grid.time < time)
  {
    solution.update();
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
  std::cout << "Execution completed in " << float(duration.count())*1e-9 << " s\n";
  std::cout << "(" << grid.iter << " iterations at "
            << float(duration.count())*1e-9/grid.iter << " s per iteration)\n";
  grid.visualize("final");
}
