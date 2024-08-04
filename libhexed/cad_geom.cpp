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

Geom::Geom(std::string file_name) {
  std::string ext = file_extension(file_name);
  HEXED_ASSERT(ext == "igs" || ext == "iges", "can only read IGES files");
  Iges_parser parser(file_name);
  auto dir = parser.section(Iges_parser::directory);
  for (auto& entry : dir) {
    Int ent_num = parser.read_int(entry[0]);
    if (ent_num == 100) {
      auto& par = parser.entry(Iges_parser::parameter, parser.read_int(entry[1]));
      std::vector<double> values;
      for (std::string s : par) values.push_back(parser.read_float(s));
      Circular_arc arc {
        {values[2], values[3], values[1]},
        std::sqrt(.5*(values[4]*values[4] + values[5]*values[5]
                      + values[6]*values[6] + values[7]*values[7])),
        std::atan2(values[5], values[4]),
        std::atan2(values[7], values[6]),
      };
      if (arc.start_angle < arc.end_angle + 1e-10) arc.start_angle += 2*constants::pi;
      arc.trans_mat = read_trans_mat(parser, parser.read_int(entry[6]));
      _curves.emplace_back(new Circular_arc {arc});
    } else if (ent_num == 124) {
    } else {
      std::cout << ent_num << "\n";
    }
  }
  Int unit_flag = parser.read_int(parser.section(Iges_parser::global)[0][13]);
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
  for (auto& curve : _curves) curve->scale = units[unit_flag - 1];
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
