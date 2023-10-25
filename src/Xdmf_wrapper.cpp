#include <Row_index.hpp>
#include <math.hpp>
#if HEXED_USE_XDMF
#include <XdmfDomain.hpp>
#include <XdmfHDF5Writer.hpp>
#include <XdmfWriter.hpp>
#include <Xdmf_wrapper.hpp>

namespace hexed
{

Xdmf_wrapper::Xdmf_wrapper(int n_dim_geom, int n_dim_topo, int row_size, std::string file_name, const Output_data& data, double time) :
  _topo{XdmfTopology::New()},
  _geom{XdmfGeometry::New()},
  _n_dim_geom{n_dim_geom},
  _n_dim_topo{n_dim_topo},
  _row_size{row_size},
  _file_name{file_name},
  _time{time},
  _n_point{math::pow(_row_size, _n_dim_topo)},
  _n_var{data.n_var(n_dim_geom)},
  _i_block{0},
  _node_inds(math::pow(2, n_dim_topo), n_dim_topo)
{
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    _attrs.push_back(XdmfAttribute::New());
    _attrs.back()->setName(data.variable_name(_n_dim_geom, i_var));
    _attrs.back()->setCenter(XdmfAttributeCenter::Node());
    _attrs.back()->setType(XdmfAttributeType::Scalar());
  }
  if (_n_dim_topo == 1) {
    _topo->setType(XdmfTopologyType::Polyline(2)); // a polyline with 2 nodes is equivalent to a line segment
    _node_inds << 0, 1;
  } else if (_n_dim_topo == 2) {
    _topo->setType(XdmfTopologyType::Quadrilateral());
    _node_inds << // XDMF uses an arbitrary node ordering, not a simple row-major order :(
      0, 0,
      1, 0,
      1, 1,
      0, 1;
  } else if (_n_dim_topo == 3) {
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
  } else HEXED_ASSERT(false, "invalid topological dimensionality");
  if      (_n_dim_geom == 2) _geom->setType(XdmfGeometryType::XY());
  else if (_n_dim_geom == 3) _geom->setType(XdmfGeometryType::XYZ());
  else HEXED_ASSERT(false, "invalid geometric dimensionality");
}

void Xdmf_wrapper::write_block(double* pos, double* vars)
{
  for (int i_elem = 0; i_elem < math::pow(_row_size - 1, _n_dim_topo); ++i_elem) {
    for (int i_vert = 0; i_vert < math::pow(2, _n_dim_topo); ++i_vert) {
      int i_node = _i_block*_n_point;
      for (int i_dim = 0; i_dim < _n_dim_topo; ++i_dim) {
        int row = (i_elem/Row_index(_n_dim_topo, _row_size - 1, i_dim).stride)%(_row_size - 1)
                  + _node_inds(i_vert, i_dim);
        i_node += row*Row_index(_n_dim_topo, _row_size, i_dim).stride;
      }
      _topo->pushBack(i_node);
    }
  }
  for (int i_point = 0; i_point < _n_point; ++i_point) {
    for (int i_dim = 0; i_dim < _n_dim_geom; ++i_dim) {
      _geom->pushBack(pos[i_dim*_n_point + i_point]);
    }
  }
  for (int i_var = 0; i_var < _n_var; ++i_var) {
    for (int i_point = 0; i_point < _n_point; ++i_point) {
      _attrs[i_var]->pushBack(vars[i_var*_n_point + i_point]);
    }
  }
  ++_i_block;
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
