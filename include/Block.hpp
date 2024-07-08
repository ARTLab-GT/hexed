#ifndef HEXED_BLOCK_HPP_
#define HEXED_BLOCK_HPP_

#include <memory>
#include "math.hpp"
#include "reciprocal.hpp"
#include "Basis.hpp"
#include "Array.hpp"
#include "Sequence.hpp"
#include "Kernel_connection.hpp"

namespace hexed::next {

class Block : public Mortal {
public:
  inline Block(int n_dim, int row_size) : _n_dim{n_dim}, _row_size{row_size} {}
  inline int n_dim() const {return _n_dim;}
  inline int row_size() const {return _row_size;}
  Mat<3> point(std::vector<int> node_coords) const;
  virtual Array<double> points() const;
  static void visualize(std::string format, std::string file_name, const std::vector<Block*>&, double time = 0.);
protected:
  virtual Mat<3> _point(std::vector<int> node_coords) const = 0;
private:
  int _n_dim;
  int _row_size;
};

class Edge;
class Element_shape;

class Vertex : public Block {
public:
  Vertex(Mat<3> pos, int row_size);
  void eat(Vertex& other);
  void glue(Element_shape& to, std::vector<double> coords);
  inline void pair(mutual::Base<Edge, Vertex>& ptr) {_edges.add(ptr);}
  inline void pair(mutual::Base<Element_shape, Vertex>& ptr) {_elems.add(ptr);}
  inline bool alive() {return _alive;}
  Mat<3> pos;
private:
  Mat<3> _point(std::vector<int>) const override;
  Reciprocal_list<Vertex, Edge> _edges;
  Reciprocal_list<Vertex, Element_shape> _elems;
  bool _alive;
  int _mass;
  Reciprocal_ptr<Vertex, Element_shape> _glued_to;
  std::vector<double> _glued_coords;
};

class Boundary_interior : public Block {
public:
  inline const Basis& basis() const {return *_basis;}
  Boundary_interior(int n_dim, const Basis& basis);
  virtual void reset() = 0;
  inline Array<double> interior() {return _interior;};
  inline bool alive() {return true;}
private:
  const Basis* _basis;
protected:
  Array<double> _interior;
};

class Edge : public Boundary_interior {
public:
  Edge(Vertex& vertex0, Vertex& vertex1, const Basis&);
  void reset() override;
  static const int no;
  void glue(Edge& other, int half = no);
  void unglue() {_glued_to.set();}
  bool glued() const {return _glued_to;}
private:
  Mat<3> _point(std::vector<int>) const override;
  std::array<Reciprocal_ptr<Edge, Vertex>, 2> _verts;
  Mortal_ptr<Edge> _glued_to;
  int _half;
};

class Face : public Boundary_interior {
public:
  Face(std::array<Vertex*, 4>, const Basis&);
  inline Edge& edge(int i) {return _edges[i];}
  void reset() override;
private:
  Mat<3> _point(std::vector<int>) const override;
  std::vector<Edge> _edges;
};

class Element_shape : public Block {
  friend class Mesh_blocks;
  friend void Vertex::glue(Element_shape&, std::vector<double>);
public:
  Element_shape(int nd, const Basis&);
  inline Vertex& vertex(int i_vert) {return *_verts[i_vert];}
  inline const Basis& basis() const {return *_basis;}
  void connect(Element_shape& other, Connection_direction);
  void connect(std::vector<Element_shape*> others, Connection_direction);
private:
  Mat<3> _vertex_point(std::vector<int>) const;
  Mat<3> _point(std::vector<int>) const override;
  const Basis* _basis;
  std::vector<Reciprocal_ptr<Element_shape, Vertex>> _verts;
  int _i_bf;
  Mortal_ptr<Boundary_interior> _bf;
  Mortal_ptr<Face> _sf;
  Reciprocal_list<Element_shape, Vertex> _glued_verts;
};

class Mesh_blocks {
public:
  Mesh_blocks(int n_dim, const Basis&);
  Sequence<Vertex&> interior_verts();
  Sequence<Vertex&> boundary_verts();
  Sequence<Edge&> edges_2d();
  Sequence<Face&> faces_3d();
  std::unique_ptr<Element_shape> create_element(Mat<3> pos, double size, int boundary_face = no_face);
  static const int no_face;
  const int n_dim;
  const Basis& basis;
private:
  std::vector<Vertex> _interior_verts;
  std::vector<Vertex> _boundary_verts;
  std::vector<Edge> _edges_2d;
  std::vector<Face> _faces_3d;
};

}
#endif
