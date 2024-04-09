#include <mesh_objects.hpp>

namespace hexed
{

Element_new::Element_new(Storage_params params, bool def, int ref_level, Array<int> pos_ind, double root_sz, Array<double> og, const Basis& b)
: _def{def},
  _params{params},
  _ref_level{ref_level},
  _pos_ind{pos_ind},
  _root_sz{root_sz},
  _origin{og},
  _basis{b},
  _vertices(std::vector<int>(_params.n_dim, 2), [this](int i) {
    Mat<> pos = nominal_position().vector();
    for (int i_dim = 0; i_dim < _params.n_dim; ++i_dim) {
      pos(i_dim) += nominal_size()*((i/math::pow(2, _params.n_dim - 1 - i_dim))%2);
    }
    return Vertex::Transferable_ptr(pos, deformed());
  }),
  //_faces({}),
  _data({}),
  tree(this)
{}

}
