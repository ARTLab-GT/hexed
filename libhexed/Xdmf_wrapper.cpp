#include <Row_index.hpp>
#include <math.hpp>
#if HEXED_USE_XDMF
#include <XdmfDomain.hpp>
#include <XdmfHDF5Writer.hpp>
#include <XdmfWriter.hpp>
#include <Xdmf_wrapper.hpp>

namespace hexed
{

Xdmf_wrapper::Xdmf_wrapper(int n_dim_geom, int n_dim_topo, std::string file_name, std::vector<std::string> var_names, double time, elem_type elem_t) :
  _topo{XdmfTopology::New()},
  _geom{XdmfGeometry::New()},
  _n_dim_geom{n_dim_geom},
  _n_dim_topo{n_dim_topo},
  _file_name{file_name},
  _time{time},
  _n_var{int(var_names.size())},
  _n_verts{0},
  _node_inds(math::pow(2, n_dim_topo), n_dim_topo),
  _elem_t{elem_t}
{
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    _attrs.push_back(XdmfAttribute::New());
    _attrs.back()->setName(var_names[i_var]);
    _attrs.back()->setCenter(XdmfAttributeCenter::Node());
    _attrs.back()->setType(XdmfAttributeType::Scalar());
  }
  HEXED_ASSERT(_n_dim_topo > 0 && _n_dim_topo <= 3, "invalid topological dimensionality");
  if (_elem_t == block) {
    if (_n_dim_topo == 1) {
      _topo->setType(XdmfTopologyType::Polyline(2)); // a polyline with 2 nodes is equivalent to a line segment
      _node_inds << 0, 1;
      _permutation.assign({0, 1});
    } else if (_n_dim_topo == 2) {
      _topo->setType(XdmfTopologyType::Quadrilateral());
      _node_inds << // XDMF uses an arbitrary node ordering, not a simple row-major order :(
        0, 0,
        1, 0,
        1, 1,
        0, 1;
      _permutation.assign({0, 2, 3, 1});
    } else {
      _topo->setType(XdmfTopologyType::Hexahedron());
      _node_inds <<
        0, 0, 0,
        1, 0, 0,
        1, 1, 0,
        0, 1, 0,
        0, 0, 1,
        1, 0, 1,
        1, 1, 1,
        0, 1, 1;
      _permutation.assign({0, 4, 6, 2, 1, 5, 7, 3});
    }
  } else {
    if      (n_dim_topo == 1) _topo->setType(XdmfTopologyType::Polyline(2));
    else if (n_dim_topo == 2) _topo->setType(XdmfTopologyType::Triangle());
    else                      _topo->setType(XdmfTopologyType::Tetrahedron());
    _permutation.resize(_n_dim_topo + 1);
    for (int i_vert = 0; i_vert < _n_dim_topo + 1; ++i_vert) _permutation[i_vert] = i_vert;
  }
  if      (_n_dim_geom == 2) _geom->setType(XdmfGeometryType::XY());
  else if (_n_dim_geom == 3) _geom->setType(XdmfGeometryType::XYZ());
  else HEXED_ASSERT(false, "invalid geometric dimensionality");
}

void Xdmf_wrapper::write_block(Array<double> pos, Array<double> vars)
{
  HEXED_ASSERT(pos.order() == _n_dim_topo + 1, "input arrays have wrong order");
  if (_n_var) {
    HEXED_ASSERT(vars.shape()[0] == _n_var, "`vars` has wrong number of rows");
    HEXED_ASSERT(pos(0).same_shape(vars(0)), "`pos` and `vars` must have compatible shape");
  }
  int row_size = pos.shape()[1];
  int n_point = math::pow(row_size, _n_dim_topo);
  for (int i_elem = 0; i_elem < math::pow(row_size - 1, _n_dim_topo); ++i_elem) {
    for (int i_vert = 0; i_vert < math::pow(2, _n_dim_topo); ++i_vert) {
      int i_node = _n_verts;
      for (int i_dim = 0; i_dim < _n_dim_topo; ++i_dim) {
        int row = (i_elem/Row_index(_n_dim_topo, row_size - 1, i_dim).stride)%(row_size - 1)
                  + _node_inds(i_vert, i_dim);
        i_node += row*Row_index(_n_dim_topo, row_size, i_dim).stride;
      }
      _topo->pushBack(i_node);
    }
  }
  for (int i_point = 0; i_point < n_point; ++i_point) {
    for (int i_dim = 0; i_dim < _n_dim_geom; ++i_dim) {
      _geom->pushBack(pos(i_dim)[i_point]);
    }
  }
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    for (int i_point = 0; i_point < n_point; ++i_point) {
      _attrs[i_var]->pushBack(vars(i_var)[i_point]);
    }
  }
  _n_verts += n_point;
}

void Xdmf_wrapper::write_unstruct(Array<int> elements, Array<double> pos, Array<double> vars)
{
  int n_elem_vert = _permutation.size();
  HEXED_ASSERT(elements.order() == 2, "`elements` must be 2D");
  HEXED_ASSERT(elements.shape()[1] == n_elem_vert, "`elements` has wrong number of columns (vertices per element)");
  HEXED_ASSERT(pos.order() == 2, "`pos` must be 2D");
  HEXED_ASSERT(pos.shape()[0] == _n_dim_geom, "`pos` must have `n_dim_geom` rows");
  int n_vert = pos.shape()[1];
  HEXED_ASSERT(vars.order() == 2, "`vars` must be 2D");
  HEXED_ASSERT(vars.shape()[0] == _n_var, "`vars` must have `n_var` rows");
  HEXED_ASSERT(vars.shape()[1] == n_vert, "`vars` and `pos` must have the same number of columns (number of vertices)");
  for (int i_elem = 0; i_elem < elements.shape()[0]; ++i_elem) {
    for (int i_vert = 0; i_vert < n_elem_vert; ++i_vert) _topo->pushBack(elements(i_elem)[_permutation[i_vert]] + _n_verts);
  }
  for (int i_vert = 0; i_vert < n_vert; ++i_vert) {
    for (int i_dim = 0; i_dim < _n_dim_geom; ++i_dim) _geom->pushBack(pos(i_dim)[i_vert]);
  }
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    for (int i_vert = 0; i_vert < n_vert; ++i_vert) _attrs[i_var]->pushBack(vars(i_var)[i_vert]);
  }
  _n_verts += pos.shape()[1];
}

Xdmf_wrapper::~Xdmf_wrapper()
{
  auto domain = XdmfDomain::New();
  auto grid = XdmfUnstructuredGrid::New();
  auto hdf5_writer = XdmfHDF5Writer::New(_file_name + ".h5");
  _topo->accept(hdf5_writer);
  _geom->accept(hdf5_writer);
  grid->setTopology(_topo);
  grid->setGeometry(_geom);
  for (auto& attr : _attrs) {
    attr->accept(hdf5_writer);
    grid->insert(attr);
  }
  grid->setTime(XdmfTime::New(_time));
  domain->insert(grid);
  domain->accept(XdmfWriter::New(_file_name + ".xmf"));
}

}
#endif
