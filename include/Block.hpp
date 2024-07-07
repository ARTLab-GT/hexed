#ifndef HEXED_BLOCK_HPP_
#define HEXED_BLOCK_HPP_

#include <memory>
#include "math.hpp"
#include "reciprocal.hpp"
#include "Basis.hpp"
#include "Array.hpp"
#include "Sequence.hpp"
#include "Kernel_connection.hpp"

namespace hexed::next
{

class Block : public Mortal
{
  protected:
  virtual Mat<3> _point(std::vector<int> node_coords) const = 0;
  public:
  const int n_dim;
  const int row_size;
  inline Block(int n_dim, int row_size) : n_dim{n_dim}, row_size{row_size} {}
  Mat<3> point(std::vector<int> node_coords) const;
  virtual Array<double> points() const;
  static void visualize(std::string format, std::string file_name, const std::vector<Block*>&, double time = 0.);
};

class Edge;
class Mesh_element;

class Vertex : public Block
{
  Reciprocal_list<Vertex, Edge> _edges;
  Reciprocal_list<Vertex, Mesh_element> _elems;
  bool _alive;
  int _mass;
  Reciprocal_ptr<Vertex, Mesh_element> _glued_to;
  std::vector<double> _glued_coords;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Mat<3> pos;
  Vertex(Mat<3> pos, int row_size);
  void eat(Vertex& other);
  void glue(Mesh_element& to, std::vector<double> coords);
  inline void pair(mutual::Base<Edge, Vertex>& ptr) {_edges.add(ptr);}
  inline void pair(mutual::Base<Mesh_element, Vertex>& ptr) {_elems.add(ptr);}
  inline bool alive() {return _alive;}
};

class Boundary_interior : public Block
{
  protected:
  Array<double> _interior;
  public:
  const Basis& basis;
  Boundary_interior(int n_dim, const Basis& basis);
  virtual void reset() = 0;
  inline Array<double> interior() {return _interior;};
  inline bool alive() {return true;}
};

class Edge : public Boundary_interior
{
  std::array<Reciprocal_ptr<Edge, Vertex>, 2> _verts;
  Mortal_ptr<Edge> _glued_to;
  int _half;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Edge(Vertex& vertex0, Vertex& vertex1, const Basis&);
  void reset() override;
  static const int no;
  void glue(Edge& other, int half = no);
  void unglue() {_glued_to.set();}
  bool glued() const {return _glued_to;}
};

class Surface_face : public Boundary_interior
{
  std::vector<Edge> _edges;
  Mat<3> _point(std::vector<int>) const override;
  public:
  Surface_face(std::array<Vertex*, 4>, const Basis&);
  inline Edge& edge(int i) {return _edges[i];}
  void reset() override;
};

class Mesh_element : public Block
{
  friend class Mesh_blocks;
  friend void Vertex::glue(Mesh_element&, std::vector<double>);
  std::vector<Reciprocal_ptr<Mesh_element, Vertex>> _verts;
  int _i_bf;
  Mortal_ptr<Boundary_interior> _bf;
  Mortal_ptr<Surface_face> _sf;
  Reciprocal_list<Mesh_element, Vertex> _glued_verts;
  Mat<3> _vertex_point(std::vector<int>) const;
  Mat<3> _point(std::vector<int>) const override;
  Mesh_element(int nd, const Basis&);
  public:
  const Basis& basis;
  inline Vertex& vertex(int i_vert) {return *_verts[i_vert];}
  void connect(Mesh_element& other, Connection_direction);
  void connect(std::vector<Mesh_element*> others, Connection_direction);
};

class Mesh_blocks
{
  std::vector<std::unique_ptr<Vertex>> _interior_verts;
  std::vector<std::unique_ptr<Vertex>> _boundary_verts;
  std::vector<std::unique_ptr<Edge>> _edges_2d;
  std::vector<std::unique_ptr<Surface_face>> _faces_3d;
  public:
  static const int no_face;
  const int n_dim;
  const Basis& basis;
  Mesh_blocks(int n_dim, const Basis&);
  Sequence<Vertex&> interior_verts();
  Sequence<Vertex&> boundary_verts();
  Sequence<Edge&> edges_2d();
  Sequence<Surface_face&> faces_3d();
  std::unique_ptr<Mesh_element> create_element(Mat<3> pos, double size, int boundary_face = no_face);
};

}
#endif
