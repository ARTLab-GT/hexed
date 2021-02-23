#ifndef CARTDG_VERTEX_HPP_
#define CARTDG_VERTEX_HPP_

#include <array>
#include <vector>

namespace cartdg
{

class Deformed_grid;

class Vertex
{
  public:
  std::array<double, 3> pos {0, 0, 0};
  std::vector<int> id_refs;
  std::vector<int> neighbor_ids;
  int mass = 1;
  int id;
  Deformed_grid* parent_grid = nullptr;
  Vertex(int id_arg);
};

}

#endif