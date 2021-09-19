#ifndef CARTDG_VERTEX_HPP_
#define CARTDG_VERTEX_HPP_

#include <array>
#include <vector>
#include <memory>
#include <set>
#include <unordered_set>

namespace cartdg
{

class Deformed_grid;

class Vertex
{
  // not thread safe!
  public:
  class Transferable_ptr;
  class Non_transferable_ptr;
  std::array<double, 3> pos {0, 0, 0};
  bool mobile = false;

  ~Vertex();
  // if we have a reason to copy/move vertices, we can implement these
  Vertex(const Vertex&) = delete;
  Vertex(Vertex&&) = delete;
  Vertex& operator=(const Vertex&) = delete;
  Vertex& operator=(Vertex&&) = delete;
  int mass();
  void eat(Vertex& other);
  void calc_relax();
  void apply_relax();

  static void connect(Vertex&, Vertex&);

  private:
  int m;
  std::array<double, 3> relax {0, 0, 0};
  std::unordered_set<Transferable_ptr*> trbl_ptrs;
  std::unordered_set<Non_transferable_ptr*> nont_ptrs;
  std::set<Vertex*> neighbors; // use set because have to iterate through
  Vertex(std::array<double, 3> pos);
};

class Vertex::Transferable_ptr
{
  std::shared_ptr<Vertex> ptr;

  public:
  Transferable_ptr(std::array<double, 3> pos);
  Transferable_ptr(const Transferable_ptr&);
  ~Transferable_ptr();
  Transferable_ptr& operator=(const Transferable_ptr&);

  Vertex* operator->();
  Vertex& operator*();
};

class Vertex::Non_transferable_ptr
{
  Vertex* ptr;

  public:
  Non_transferable_ptr(Vertex&);
  Non_transferable_ptr(const Non_transferable_ptr&);
  ~Non_transferable_ptr();
  Non_transferable_ptr& operator=(const Non_transferable_ptr&);

  operator bool();
  Vertex* operator->();
  Vertex& operator*();
  void nullify();
};

}

#endif
