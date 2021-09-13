#ifndef CARTDG_VERTEX_HPP_
#define CARTDG_VERTEX_HPP_

#include <array>
#include <vector>
#include <memory>
#include <unordered_set>

namespace cartdg
{

class Deformed_grid;

class Vertex
{
  public:
  class Non_transferable_ptr;
  class Transferable_ptr;
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

  private:
  int m;
  std::unordered_set<Transferable_ptr*> trbl_ptrs;
  Vertex(std::array<double, 3> pos);
};

class Vertex::Non_transferable_ptr
{
  public:
  Non_transferable_ptr(Vertex& target);
  Non_transferable_ptr(const Non_transferable_ptr& other);
  Non_transferable_ptr(Non_transferable_ptr&& other);
  ~Non_transferable_ptr();
  Non_transferable_ptr& operator=(const Non_transferable_ptr& other);
  Non_transferable_ptr& operator=(Non_transferable_ptr&& other);

  bool is_assigned();
  Vertex* operator->();
};

class Vertex::Transferable_ptr
{
  std::shared_ptr<Vertex> ptr;

  public:
  Transferable_ptr(std::array<double, 3> pos);
  Transferable_ptr(const Transferable_ptr&);
  Transferable_ptr(Transferable_ptr&&);
  ~Transferable_ptr();
  Transferable_ptr& operator=(const Transferable_ptr&);
  Transferable_ptr& operator=(Transferable_ptr&&);

  Vertex* operator->();
  Vertex& operator*();
};

}

#endif
