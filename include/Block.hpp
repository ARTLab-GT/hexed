#ifndef HEXED_BLOCK_HPP_
#define HEXED_BLOCK_HPP_

#include <memory>
#include "math.hpp"
#include "Mutual_ptr.hpp"
#include "Basis.hpp"
#include "Array.hpp"

namespace hexed::next
{

class Block
{
  protected:
  virtual Mat<3> _point(std::vector<int> node_coords) const = 0;
  public:
  const int n_dim;
  const int row_size;
  inline Block(int n_dim, int row_size) : n_dim{n_dim}, row_size{row_size} {}
  Mat<3> point(std::vector<int> node_coords) const;
};

class Edge;

class Vertex : public Block
{
  std::vector<Mutual_ptr<Vertex, Edge>> _edges;
  inline Mat<3> _point(std::vector<int>) const override {return pos;}
  bool _alive;
  int _mass;
  public:
  Mat<3> pos;
  inline Vertex(Mat<3> pos, int row_size) : Block{0, row_size}, pos{pos}, _alive{true}, _mass{1} {}
  void eat(Vertex& other);
  void average(std::vector<Vertex*>);
  void pair(Mutual_ptr<Edge, Vertex>& ptr);
  void purge();
  inline bool alive() {return _alive;}
};

class Edge : public Block
{
  int _rs;
  std::array<Mutual_ptr<Edge, Vertex>, 2> _verts;
  Array<double> _interior;
  std::shared_ptr<Basis> _basis;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis>);
  inline Array<double> interior() {return _interior;};
  void reset();
  void glue(Edge& other, double start, double stop);
};

}
#endif
