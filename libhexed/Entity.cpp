#include <Entity.hpp>

namespace hexed
{

Entity::Entity(int nd, int md, const std::shared_ptr<Basis> b)
: n_dim{nd}, my_dim{nd}, basis{b}
{}

Cartesian::Cartesian(int nd, int md, const std::shared_ptr<Basis> b, Array<int> nom_pos, double nom_sz)
: Entity{nd, md, b}, _nom_pos(nom_pos.copy()), _nom_sz{nom_sz}
{}

Array<int> Cartesian::nominal_position() const {return _nom_pos.copy();}
double Cartesian::nominal_size() const {return _nom_sz;}

const Array<double> Cartesian::points() const
{
  std::vector<int> sz{n_dim};
  int rs = basis->row_size;
  auto nodes = basis->nodes();
  for (int i_dim = 0; i_dim < my_dim; ++i_dim) sz.push_back(rs);
  Array<double> pts(sz);
  for (int i_pt = 0; i_pt < pts(0).size(); ++i_pt) {
    for (int i_dim = 0; i_dim < my_dim; ++i_dim) {
      int i_row = (i_pt/pts.stride(1 + i_dim))%rs;
      pts(i_dim)[i_pt] = _nom_sz*(_nom_pos[i_dim] + nodes(i_row));
    }
  }
  return pts;
}

}
