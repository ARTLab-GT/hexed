#include <hexed/cad_geom.hpp>
#include <hexed/Iges_parser.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/utils.hpp>
#include <hexed/constants.hpp>

namespace hexed::cad_geom {

Mat<3> Circular_arc::temp_point(Mat<1> params) const {
  double angle = start_angle + params(0)*(start_angle - end_angle);
  return center + radius*Mat<3>{std::cos(angle), std::sin(angle), 0.};
}

Revolution_surface::Revolution_surface(Entity<1>* g, Line_segment ax, double sa, double ea)
: generatrix{g}, axis{ax}, start_angle{sa}, end_angle{ea}
{}

Mat<3> Revolution_surface::temp_point(Mat<2> params) const {
  Mat<3> unrotated = generatrix->point(params(Eigen::seqN(0, 1))) - axis.endpoints(all, 0);
  double angle = start_angle + params(1)*(start_angle - end_angle);
  Mat<3> ax_vec = (axis.endpoints(all, 1) - axis.endpoints(all, 0)).normalized();
  Mat<3> axial_component = unrotated.dot(ax_vec)*ax_vec;
  Mat<3> radial_component = unrotated - axial_component;
  return axial_component + std::cos(angle)*radial_component + std::sin(angle)*ax_vec.cross(radial_component);
}

Trans_mat read_trans_mat(const Iges_parser& parser, Int line) {
  Trans_mat tm;
  if (line > 0) {
    const auto& dir_entry = parser.entry(Iges_parser::directory, line);
    HEXED_ASSERT(parser.read_int(dir_entry[6]) == 0, "chaining transformation matrices is not yet implemented");
    auto& par_entry = parser.entry(Iges_parser::parameter, parser.read_int(dir_entry[1]));
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      for (int j_dim = 0; j_dim < 3; ++j_dim) {
        tm.transform(i_dim, j_dim) = parser.read_float(par_entry[4*i_dim + j_dim + 1]);
      }
      tm.translate(i_dim) = parser.read_float(par_entry[4*i_dim + 3 + 1]);
    }
  }
  return tm;
}

template <typename T>
class Read_entity {
  public:
  Read_entity(const Iges_parser& parser, Int line)
  : Read_entity{parser, parser.entry(Iges_parser::directory, line)}
  {}
  Read_entity(const Iges_parser& parser, const std::vector<std::string>& dir)
  : _parser{parser}
  , _dir{dir}
  , _par{_parser.entry(Iges_parser::parameter, _parser.read_int(_dir[1]))}
  , _ent_num{_parser.read_int(_dir[0])}
  {}

  void read_circular_arc() {
    if (_ptr || _ent_num != 100) return;
    std::vector<double> values;
    for (std::string s : _par) values.push_back(_parser.read_float(s));
    auto arc = new Circular_arc {
      {values[2], values[3], values[1]},
      std::sqrt(.5*(values[4]*values[4] + values[5]*values[5]
                    + values[6]*values[6] + values[7]*values[7])),
      std::atan2(values[5], values[4]),
      std::atan2(values[7], values[6]),
    };
    if (arc->start_angle < arc->end_angle + 1e-10) arc->start_angle += 2*constants::pi;
    _ptr.reset(arc);
  }

  void read_line_segment() {
    if (_ptr || _ent_num != 110) return;
    Mat<3, 2> endpts;
    for (int i = 0; i < 6; ++i) endpts(i) = _parser.read_float(_par[1 + i]);
    _ptr.reset(new Line_segment(endpts));
  }

  void read_curve() {
    read_circular_arc();
    read_line_segment();
  }

  std::unique_ptr<T> get() {
    if (_ptr) {
      _ptr->trans_mat = read_trans_mat(_parser, _parser.read_int(_dir[6]));
      Int unit_flag = _parser.read_int(_parser.section(Iges_parser::global)[0][13]);
      HEXED_ASSERT(unit_flag > 0 && unit_flag <= 11, "invalid unit flag");
      double units [] {
        constants::inch,
        1e-3*constants::meter,
        1.,
        constants::foot,
        constants::mile,
        constants::meter,
        1e3*constants::meter,
        1e-3*constants::inch,
        1e-6*constants::meter,
        1e-2*constants::meter,
        1e-6*constants::inch,
      };
      _ptr->scale = units[unit_flag - 1];
    }
    return std::unique_ptr<T>{_ptr.release()};
  }

  Int entity_number() const {return _ent_num;}

  private:
  const Iges_parser& _parser;
  const std::vector<std::string>& _dir;
  const std::vector<std::string>& _par;
  Int _ent_num;
  std::unique_ptr<T> _ptr;
};

Geom::Geom(std::string file_name) {
  std::string ext = file_extension(file_name);
  HEXED_ASSERT(ext == "igs" || ext == "iges", "can only read IGES files");
  Iges_parser parser(file_name);
  auto dir = parser.section(Iges_parser::directory);
  for (auto& entry : dir) {
    Read_entity<Entity<1>> read(parser, entry);
    read.read_curve();
    std::unique_ptr<Entity<1>> ptr = read.get();
    if (ptr) _curves.emplace_back(ptr.release());
    else std::cout << read.entity_number() << "\n";
  }
}

void Geom::visualize(std::string file_name) const {
  auto vis = Visualizer::create("default", 3, 1, "edges", {}, 0., Visualizer::block);
  int n = 101;
  for (auto& c : _curves) {
    Array<double> discrete({3, n});
    for (int i = 0; i < n; ++i) {
      Mat<3> p = c->point(Mat<1>{i/(n - 1.)});
      for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)[i] = p(i_dim);
    }
    vis->write_block(discrete, Array<double>({0, n}));
  }
}

}
