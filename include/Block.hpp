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
  Mutual_ptr<Edge, Edge> _glued_to;
  int _half;
  std::vector<Mutual_ptr<Edge, Edge>> _glued;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis>);
  inline Array<double> interior() {return _interior;};
  void reset();
  static const int no;
  void glue(Edge& other, int half = no);
  void unglue();
  bool glued() const;
};

class Surface_face : public Block
{
  int _rs;
  std::vector<Edge> _edges;
  Array<double> _interior;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Surface_face(std::array<Vertex*, 4>, std::shared_ptr<Basis>);
  inline Array<double> interior() {return _interior;};
  inline Edge& edge(int i) {return _edges[i];}
  void reset();
};

}
#endif
