#include <iostream>
#include <chrono>

#include <Solution.hpp>

int main()
{
  // Resolution parameters
  const int rank = 6;
  const int ref_level = 3;

  // Solution setup
  cartdg::Solution solution (4, 2, rank, 1.);
  solution.add_block_grid(ref_level);
  cartdg::Grid& grid = solution.get_grid(0);
  int n_div = 1; for (int i = 0; i < ref_level; ++i) n_div *= 2;
  std::vector<int> periods {n_div, n_div};
  grid.auto_connect(periods);
  cartdg::Isentropic_vortex vortex (std::vector<double> {10., 0., 1.225, 101325/0.4});
  solution.initialize(vortex);

  // Let's go!
  solution.visualize("initial");
  auto start = std::chrono::high_resolution_clock::now();
  double time = 0;
  for (int i = 0; i < 20; ++i)
  {
    time += 0.001;
    while (grid.time < time)
    {
      solution.update();
    }
    solution.visualize("simulation");
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
  solution.visualize("final");
}
