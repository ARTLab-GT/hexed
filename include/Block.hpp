#include <memory>
#include "math.hpp"
#include "Mutual_ptr.hpp"
#include "Basis.hpp"

namespace hexed::next
{

class Block
{
  public:
  const int n_dim;
  const int row_size;
  inline Block(int n_dim, int row_size) : n_dim{n_dim}, row_size{row_size} {}
  virtual Mat<3> point(std::vector<int> node_coords) const = 0;
};

class Vertex : public Block
{
  std::vector<Mutual_ptr<Vertex, Edge>> _edges;
  public:
  Mat<3> pos;
  inline Vertex(Mat<3> pos, int row_size) : Block(0, row_size), pos{pos} {}
  inline Mat<3> point(std::vector<int>) const override {return pos;}
  void eat(Vertex& other);
  void average(std::vector<Vertex*>);
};

class Edge : public Block
{
  std::array<Mutual_ptr<Edge, Vertex>, 2> _vertices;
  Array<double> _interior;
  std::shared_ptr<Basis> basis;
  public:
  Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis>);
  Mat<3> point(std::vector<int>) const override;
  Array<double> interior();
  void glue(Edge& other, double start, double stop);
};

}
