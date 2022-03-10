#include <iostream>
#include <chrono>

#include <Solution.hpp>

int main()
{
  // Resolution parameters
  const int row_size = 6;
  const int ref_level = 3;

  // Solution setup
  cartdg::Solution solution (4, 2, row_size, 1.);
  solution.add_block_grid(ref_level);
  cartdg::Regular_grid& grid = solution.reg_grids[0];
  int n_div = 1; for (int i = 0; i < ref_level; ++i) n_div *= 2;
  std::vector<int> periods {n_div, n_div};
  grid.auto_connect(periods);
  cartdg::Isentropic_vortex vortex (std::vector<double> {100., 0., 1.225, 101325/0.4});
  vortex.center0 = 0.5; vortex.center1 = 0.5;
  solution.initialize(vortex);

  // Let's go!
  solution.visualize_field("demo_initial");
  auto start = std::chrono::high_resolution_clock::now();
  double time = 0;
  for (int i = 0; i < 20; ++i)
  {
    time += 1e-3;
    while (grid.time < time)
    {
      solution.update(0.9);
    }
    char buffer [100];
    snprintf(buffer, 100, "demo_time_%e", time);
    solution.visualize_field(buffer);
  }
  time = 0.2;
  while (grid.time < time)
  {
    solution.update(0.9);
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
  std::cout << "Execution completed in " << float(duration.count())*1e-9 << " s\n";
  std::cout << "(" << grid.iter << " iterations at "
            << float(duration.count())*1e-9/grid.iter << " s per iteration)\n";
  solution.visualize_field("demo_final");
}
