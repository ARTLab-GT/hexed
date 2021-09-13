#ifndef CARTDG_VERTEX_HPP_
#define CARTDG_VERTEX_HPP_

#include <array>
#include <vector>
#include <memory>

namespace cartdg
{

class Deformed_grid;

class Vertex
{
  int m = 1;
  Vertex(std::array<double, 3> pos);

  public:
  class Non_transferrable_ptr;
  class Transferrable_ptr;
  std::array<double, 3> pos {0, 0, 0};
  bool mobile = false;

  Vertex(int id_arg);
  ~Vertex();
  // FIXME: handle copy/move
  int mass();
  void eat(Vertex& other);
  void calc_relax();
  void apply_relax();

  // should not be public
  std::array<double, 3> relax {0, 0, 0};
  std::vector<int> id_refs;
  std::vector<int> neighbor_ids;
  int id;
  Deformed_grid* parent_grid = nullptr;
};

class Vertex::Non_transferrable_ptr
{
  public:
  Non_transferrable_ptr(Vertex& target);
  Non_transferrable_ptr(const Non_transferrable_ptr& other);
  Non_transferrable_ptr(Non_transferrable_ptr&& other);
  ~Non_transferrable_ptr();
  Non_transferrable_ptr& operator=(const Non_transferrable_ptr& other);
  Non_transferrable_ptr& operator=(Non_transferrable_ptr&& other);

  bool is_assigned();
  Vertex* operator->();
};

class Vertex::Transferrable_ptr
{
  std::shared_ptr<Vertex> ptr;

  public:
  Transferrable_ptr(std::array<double, 3> pos);
  Transferrable_ptr(const Transferrable_ptr& other);
  Transferrable_ptr(Transferrable_ptr&& other);
  ~Transferrable_ptr();
  Transferrable_ptr& operator=(const Transferrable_ptr& other);
  Transferrable_ptr& operator=(Transferrable_ptr&& other);

  Vertex* operator->();
  Vertex& operator*();
};

}

#endif
