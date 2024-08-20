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
{
  HEXED_ASSERT(end_angle - start_angle > 0, "end angle must be greater than start angle");
}

Mat<2> Revolution_surface::temp_nearest_params(Mat<3> p) const {
  Mat<3> unit_axis = (axis.endpoints(all, 1) - axis.endpoints(all, 0)).normalized();
  Mat<3> from_start = p - axis.endpoints(all, 0);
  Mat<3> radius = (from_start - from_start.dot(unit_axis)*unit_axis).normalized();
  double dist = huge;
  Mat<2> nearest {std::nan(""), std::nan("")};
  for (int i = 0; i < n_div + 1; ++i) {
    double param = i/double(n_div);
    Mat<3> arc_point = generatrix->point(Mat<1>{param}) - axis.endpoints(all, 0);
    Mat<3> arc_radius = (arc_point - arc_point.dot(unit_axis)*unit_axis).normalized();
    double angle = std::atan2(arc_radius.cross(radius).dot(unit_axis), arc_radius.dot(radius));
    if (math::angle_diff(angle, start_angle) > end_angle - start_angle) {
      angle = (math::angle_diff(angle, end_angle) > math::angle_diff(start_angle, angle)) ? start_angle : end_angle;
    }
    Mat<2> candidate {param, math::angle_diff(angle, start_angle)/(end_angle - start_angle)};
    double d = (temp_point(candidate) - p).norm();
    if (d < dist) {
      dist = d;
      nearest = candidate;
    }
  }
  return nearest;
}

Mat<3> Revolution_surface::temp_point(Mat<2> params) const {
  Mat<3> unrotated = generatrix->point(params(Eigen::seqN(0, 1))) - axis.endpoints(all, 0);
  double angle = start_angle + params(1)*(start_angle - end_angle);
  Mat<3> ax_vec = (axis.endpoints(all, 1) - axis.endpoints(all, 0)).normalized();
  Mat<3> axial_component = unrotated.dot(ax_vec)*ax_vec;
  Mat<3> radial_component = unrotated - axial_component;
  return axial_component + std::cos(angle)*radial_component + std::sin(angle)*radial_component.cross(ax_vec);
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
    //HEXED_ASSERT(_ptr, format_str(200, "could not read entity type `%i` as curve", int(_ent_num)));
  }

  void read_plane() {
    if (_ptr || _ent_num != 108) return;
    HEXED_ASSERT(_parser.read_int(_par[5]) == 0, "bounded planes are not implemented", assert::Not_implemented_error);
    double coefs [4];
    for (int i_coef = 0; i_coef < 4; ++i_coef) coefs[i_coef] = _parser.read_float(_par[i_coef + 1]);
    auto comp = [](double x, double y){return std::abs(x) < std::abs(y);};
    int i_dependent = std::max_element(coefs, coefs + 3, comp) - coefs;
    int vec_inds [2] {(i_dependent + 1)%3, (i_dependent + 2)%3};
    Mat<3> origin = Mat<3>::Zero();
    origin(i_dependent) = coefs[3]/coefs[i_dependent];
    Mat<3, 2> vecs = Mat<3, 2>::Zero();
    for (int i_vec = 0; i_vec < 2; ++i_vec) {
      vecs(vec_inds[i_vec], i_vec) = 1;
      vecs(i_dependent, i_vec) = -coefs[vec_inds[i_vec]]/coefs[i_dependent];
    }
    _ptr.reset(new Plane(origin, vecs*1000));
  }

  void read_revolution_surface() {
    if (_ptr || _ent_num != 120) return;
    Read_entity<Line_segment> read_axis(_parser, _parser.read_int(_par[1]));
    read_axis.read_line_segment();
    Read_entity<Entity<1>> read_generatrix(_parser, _parser.read_int(_par[2]));
    read_generatrix.read_curve();
    auto surf = new Revolution_surface {
      read_generatrix.get().release(),
      *read_axis.get(),
      _parser.read_float(_par[3]),
      _parser.read_float(_par[4]),
    };
    surf->generatrix->scale = 1.;
    surf->axis.scale = 1.;
    _ptr.reset(surf);
  }

  void read_surface() {
    read_plane();
    read_revolution_surface();
  }

  void read_disc_curve() {
    Read_entity<Entity<1>> curve(_parser, _dir);
    curve.read_curve();
    auto ptr = curve.get();
    _ptr.reset(new Composite_curve);
    if (ptr) _ptr->emplace_back(ptr.release());
    else if (curve.entity_number() == 116 || curve.entity_number() == 132) {}
    else if (curve.entity_number() == 102) {
      for (int i_curve = 0; i_curve < _parser.read_int(_par[1]); ++i_curve) {
        Read_entity<T> sub_curve(_parser, _parser.read_int(_par[2 + i_curve]));
        sub_curve.read_disc_curve();
        for (auto& c : *sub_curve._ptr) _ptr->emplace_back(c.release());
      }
    } else HEXED_THROW(format_str(200, "could not read entity type `%i` as a curve", int(curve.entity_number())));
  }

  void read_trimmed_surface() {
    if (_ptr || _ent_num != 144) return;
    Read_entity<Entity<2>> read_surf(_parser, _parser.read_int(_par[1]));
    read_surf.read_surface();
    _ptr.reset(new Trimmed_surface());
    _ptr->surface.reset(read_surf.get().release());
    if (!_ptr->surface) {
      std::cout << format_str(100, "failed to read surface from entity # %lli\n", read_surf.entity_number());
    } else _ptr->surface->scale = 1.;
    HEXED_ASSERT(_parser.read_int(_par[2]) == 1, "using outer boundary as boundary curve is not implemented",
                 assert::Not_implemented_error);
    Read_entity<T> on_surf(_parser, _parser.read_int(_par[4]));
    HEXED_ASSERT(on_surf.entity_number() == 142, "boundary must be a curve on a surface");
    int model_curve = _parser.read_int(on_surf._par[4]);
    int prev_div = -1;
    double prev_ordinate;
    if ((_parser.read_int(on_surf._par[5]) == 1 || model_curve == 0) && _parser.read_int(on_surf._par[3]) != 0) {
      HEXED_THROW("parameter-space curves are not implemented", assert::Not_implemented_error);
    } else {
      HEXED_ASSERT(model_curve != 0,
                   "at least one of parameter-space and model-space curve pointers must be defined");
      Read_entity<Composite_curve> curve(_parser, model_curve);
      curve.read_disc_curve();
      _ptr->curves.emplace_back(curve.get().release());
      double sz = 1./n_div;
      for (int i_div = 0; i_div < n_div; ++i_div) {
        Mat<2, 2> node_params;
        for (int endpoint = 0; endpoint < 2; ++endpoint) {
           node_params(all, endpoint) = _ptr->surface->nearest_params(
             _ptr->curves.back()(Mat<1>{(i_div + endpoint)*sz})
           );
        }
        if (prev_div < 0);
          prev_div = 
          prev_param = node_param(0, 1);
        }
        while (std::abs(node_params(0, 1) - prev_params(0)) >= sz
               && std::abs(node_params(0, 1) - node_params(0, 0)) > 1e-12) {
          bool sign = node_params(0, 1) > prev_params(0);
          prev_params(0) += math::sign(sign)*sz;
          prev_params(1) = node_params(1, 0) + (prev_params(0) - node_params(0, 0))
                                                *(node_params(1, 1) - node_params(1, 0))
                                                /(node_params(0, 1) - node_params(0, 0));
          node_params(all, 0) = prev_params;
        }
      }
    }
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

template <> std::unique_ptr<Composite_curve> Read_entity<Composite_curve>::get() {
  return std::unique_ptr<Composite_curve>{_ptr.release()};
}

Geom::Geom(std::string file_name) {
  std::string ext = file_extension(file_name);
  HEXED_ASSERT(ext == "igs" || ext == "iges", "can only read IGES files");
  Iges_parser parser(file_name);
  auto dir = parser.section(Iges_parser::directory);
  for (auto& entry : dir) {
    Read_entity<Trimmed_surface> read(parser, entry);
    read.read_trimmed_surface();
    std::unique_ptr<Trimmed_surface> ts(read.get());
    if (ts) _surfaces.emplace_back(ts.release());
  }
}

void Geom::visualize(std::string file_name) const {
  int n = 61;
  {
    auto vis = Visualizer::create("default", 3, 1, "edges", {}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      for (auto& c : s->curves) {
        for (auto& curve : *c) {
          Array<double> discrete({3, n});
          for (int i = 0; i < n; ++i) {
            Mat<3> p = curve->point(Mat<1>{i/(n - 1.)});
            for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)[i] = p(i_dim);
          }
          vis->write_block(discrete, Array<double>({0, n}));
        }
      }
    }
  }
  {
    auto vis = Visualizer::create("default", 3, 2, "surfaces", {}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      Array<double> discrete({3, n, n});
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          Mat<3> p = s->point(Mat<2>{i/(n - 1.), j/(n - 1.)});
          for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)(i)[j] = p(i_dim);
        }
      }
      vis->write_block(discrete, Array<double>({0, n, n}));
    }
  }
  {
    auto vis = Visualizer::create("default", 3, 3, "distance", {"distance"}, 0., Visualizer::block);
    Array<double> coords({3, n, n, n});
    Array<double> dist({1, n, n, n});
    for (int i = 0; i < math::pow(n, 3); ++i) {
      Mat<3> p;
      for (int i_dim = 0; i_dim < 3; ++i_dim) p(i_dim) = coords(i_dim)[i] = 1./n*(i/math::pow(n, 2 - i_dim)%n);
      double min_dist = huge;
      for (auto& s : _surfaces) min_dist = std::min(min_dist, (p - s->nearest_point(p)).norm());
      dist[i] = min_dist;
    }
    vis->write_block(coords, dist);
  }
}

}
