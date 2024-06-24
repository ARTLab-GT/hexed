#ifndef HEXED_BLOCK_HPP_
#define HEXED_BLOCK_HPP_

#include <memory>
#include "math.hpp"
#include "Mutual_ptr.hpp"
#include "Basis.hpp"
#include "Array.hpp"
#include "Sequence.hpp"

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
  virtual ~Block() = default;
  Mat<3> point(std::vector<int> node_coords) const;
  virtual Array<double> points() const;
  static void visualize(std::string format, std::string file_name, const std::vector<Block*>&, double time = 0.);
};

class Edge;
class Mesh_element;

class Vertex : public Block
{
  std::vector<Mutual_ptr<Vertex, Edge>> _edges;
  std::vector<Mutual_ptr<Vertex, Mesh_element>> _elems;
  inline Mat<3> _point(std::vector<int>) const override {return pos;}
  bool _alive;
  int _mass;
  public:
  Mat<3> pos;
  inline Vertex(Mat<3> pos, int row_size) : Block{0, row_size}, pos{pos}, _alive{true}, _mass{1} {}
  void eat(Vertex& other);
  void pair(Mutual_ptr<Edge, Vertex>& ptr);
  void pair(Mutual_ptr<Mesh_element, Vertex>& ptr);
  void purge();
  inline bool alive() {return _alive;}
};

class Boundary_interior : public Block
{
  protected:
  std::shared_ptr<Basis> _basis;
  Array<double> _interior;
  public:
  Boundary_interior(int n_dim, std::shared_ptr<Basis> basis);
  virtual void reset() = 0;
  inline Array<double> interior() {return _interior;};
};

class Edge : public Boundary_interior
{
  std::array<Mutual_ptr<Edge, Vertex>, 2> _verts;
  Mutual_ptr<Edge, Edge> _glued_to;
  int _half;
  std::vector<Mutual_ptr<Edge, Edge>> _glued;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Edge(Vertex& vertex0, Vertex& vertex1, std::shared_ptr<Basis>);
  void reset() override;
  static const int no;
  void glue(Edge& other, int half = no);
  void unglue();
  bool glued() const;
};

class Surface_face : public Boundary_interior
{
  std::vector<Edge> _edges;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Surface_face(std::array<Vertex*, 4>, std::shared_ptr<Basis>);
  inline Edge& edge(int i) {return _edges[i];}
  void reset() override;
};

class Mesh_element : public Block
{
  friend class Mesh_blocks;
  std::vector<Mutual_ptr<Mesh_element, Vertex>> _verts;
  Mat<3> _point(std::vector<int>) const override;
  Mesh_element(int nd, int rs);
  public:
  inline Vertex& vertex(int i_vert) {return *_verts[i_vert];}
};

class Mesh_blocks
{
  std::shared_ptr<Basis> _basis;
  std::vector<std::unique_ptr<Vertex>> _interior_verts;
  std::vector<std::unique_ptr<Vertex>> _boundary_verts;
  std::vector<std::unique_ptr<Edge>> _2d_edges;
  std::vector<std::unique_ptr<Surface_face>> _3d_faces;
  public:
  static const int no_face;
  const int n_dim;
  Mesh_blocks(int n_dim, std::shared_ptr<Basis>);
  inline const Basis& basis() const {return *_basis;}
  Sequence<Vertex&> interior_verts();
  Sequence<Vertex&> boundary_verts();
  std::unique_ptr<Mesh_element> create_element(Mat<3> pos, double size, int boundary_face = no_face);
};

}
#endif
