#include <hexed/brep.hpp>
#include <hexed/Iges_parser.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/utils.hpp>
#include <hexed/constants.hpp>

namespace hexed::brep {

Entity<1>::Nearest_params Line_segment::temp_nearest_params(
  Mat<3> p, Entity<1>::Constraint is_feasible, double max_distance
) const {
  Mat<3> diff = endpoints(all, 1) - endpoints(all, 0);
  Mat<1> params {std::max(0., std::min(1., (p - endpoints(all, 0)).dot(diff)/diff.squaredNorm()))};
  return {params, is_feasible(params)};
}

double limited_angle(double angle, double start, double end) {
  if (math::angle_diff(angle, start) > end - start) {
    return (math::angle_diff(angle, end) > math::angle_diff(start, angle)) ? start : end;
  }
  return angle;
}

Entity<1>::Nearest_params Circular_arc::temp_nearest_params(
  Mat<3> p, Entity<1>::Constraint is_feasible, double max_distance
) const {
  double angle = limited_angle(std::atan2(p(1) - center(1), p(0) - center(0)), start_angle, end_angle);
  Mat<1> params {math::angle_diff(angle, start_angle)/(end_angle - start_angle)};
  return {params, is_feasible(params)};
}

Mat<3> Circular_arc::temp_point(Mat<1> params) const {
  double angle = start_angle + params(0)*(end_angle - start_angle);
  return center + radius*Mat<3>{std::cos(angle), std::sin(angle), 0.};
}

Array<double> discretize(Entity<1>& curve, Int n_div) {
  Array<double> nodes({n_div + 1, 3});
  for (Int i_node = 0; i_node < n_div + 1; ++i_node) nodes(i_node).vector() = curve.point(Mat<1>{i_node/double(n_div)});
  return nodes;
}

Revolution_surface::Revolution_surface(Entity<1>* g, Line_segment ax, double sa, double ea)
: generatrix{g}, axis{ax}, start_angle{sa}, end_angle{ea}, _tree(discretize(*generatrix, n_div), 4)
{
  HEXED_ASSERT(end_angle - start_angle > 0, "end angle must be greater than start angle");
}

Mat<3> Revolution_surface::rotate(Mat<3> p, double angle) const {
  p -= axis.point(Mat<1>{0.});
  Mat<3> ax_vec = (axis.point(Mat<1>{1.}) - axis.point(Mat<1>{0.})).normalized();
  Mat<3> axial_component = p.dot(ax_vec)*ax_vec;
  Mat<3> radial_component = p - axial_component;
  return std::cos(angle)*radial_component + std::sin(angle)*radial_component.cross(ax_vec)
         + axial_component + axis.point(Mat<1>{0.});
};

class Revo_surf_tnp {
  public:
  struct Candidate {
    Entity<2>::Nearest_params np;
    double dist;
  };
  Revo_surf_tnp(const Revolution_surface& s, Mat<3> p, Entity<2>::Constraint is_f, double max_distance)
  : surf{s}
  , is_feasible{is_f}
  , point{p}
  , unit_axis{(surf.axis.point(Mat<1>{1.}) - surf.axis.point(Mat<1>{0.})).normalized()}
  , from_start{p - surf.axis.point(Mat<1>{0.})}
  , radius{(from_start - from_start.dot(unit_axis)*unit_axis).normalized()}
  , cand{{Mat<2>{std::nan(""), std::nan("")}, false}, max_distance}
  {}
  double best_angle(Mat<3> arc_point) {
    arc_point -= surf.axis.point(Mat<1>{0.});
    Mat<3> arc_radius = (arc_point - arc_point.dot(unit_axis)*unit_axis).normalized();
    double angle = std::atan2(arc_radius.cross(radius).dot(unit_axis), arc_radius.dot(radius));
    return limited_angle(angle, surf.start_angle, surf.end_angle);
  }
  Mat<3> best_point(Mat<3> arc_point) {
    double angle = best_angle(arc_point);
    return surf.rotate(arc_point, angle);
  }
  Candidate merge(Candidate c0, Candidate c1) {
    if (c1.dist < c0.dist && c1.np.is_feasible) return c1;
    return c0;
  }
  void find(const Tree_curve::Segment& segment) {
    #if 0
    if ((best_point(segment.center) - point).norm() - segment.radius < cand.dist) {
      if (segment.segments.size()) {
        for (auto& seg : segment.segments) find(seg);
      } else {
    #endif
        Int n_nodes = segment.nodes.shape()[0];
        for (Int i_node = 0; i_node < n_nodes; ++i_node) {
          Candidate c;
          Mat<3> node = segment.nodes(i_node).vector();
          c.np.params(0) = double(segment.nodes_start + i_node)/surf.n_div;
          double angle = best_angle(node);
          c.np.params(1) = math::angle_diff(angle, surf.start_angle)/(surf.end_angle - surf.start_angle);
          c.np.is_feasible = is_feasible(c.np.params);
          c.dist = (surf.rotate(node, angle) - point).norm();
          cand = merge(cand, c);
        }
    #if 0
      }
    }
    #endif
  }
  const Revolution_surface& surf;
  Entity<2>::Constraint is_feasible;
  Mat<3> point;
  Mat<3> unit_axis;
  Mat<3> from_start;
  Mat<3> radius;
  Candidate cand;
};

Entity<2>::Nearest_params Revolution_surface::temp_nearest_params(
  Mat<3> p, Entity<2>::Constraint is_feasible, double max_distance
) const {
  Revo_surf_tnp tnp(*this, p, is_feasible, max_distance);
  tnp.find(_tree.root());
  return tnp.cand.np;
}

Mat<3> Revolution_surface::temp_point(Mat<2> params) const {
  Mat<3> p = generatrix->point(params(Eigen::seqN(0, 1)));
  double angle = start_angle + params(1)*(start_angle - end_angle);
  return rotate(p, angle);
}

bool Trimmed_surface::inside(Mat<2> params) const {
  Int i_seg = floor(params(0)*n_div);
  if (i_seg < 0 || i_seg > n_div) return false;
  if (i_seg == n_div) i_seg = n_div - 1;
  auto& segments = parametric_segments[i_seg];
  Int n_intersections = 0;
  for (Mat<2> seg : segments) n_intersections += params(1) < seg(0) + (params(0)*n_div - i_seg)*(seg(1) - seg(0));
  return n_intersections%2;
}

Entity<2>::Nearest_params Trimmed_surface::temp_nearest_params(
  Mat<3> p, Entity<2>::Constraint is_feasible, double max_distance
) const {
  auto f = [this, is_feasible](Mat<2> params){return inside(params) && is_feasible(params);};
  Nearest_params nearest = surface->nearest_params(p, f, max_distance);
  double dist = (temp_point(nearest.params) - p).norm();
  for (auto& composite : curves) {
    for (auto& curve : *composite) {
      Mat<3> candidate = curve->nearest_point(p, max_distance);
      double d = (candidate - p).norm();
      if (d < dist || !nearest.is_feasible) {
        Nearest_params par = surface->nearest_params(candidate, is_feasible);
        // note we use `is_feasible` instead of `f` because a point on the boundary is guaranteed to be inside,
        // even if (especially if) the parametric boundary segments mark it as outside
        if (par.is_feasible) {
          nearest.params = par.params;
          nearest.is_feasible = true;
          dist = d;
        }
      }
    }
  }
  return nearest;
};

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
    _ptr.reset(new Plane(origin, vecs));
  }

  void read_revolution_surface() {
    if (_ptr || _ent_num != 120) return;
    Read_entity<Line_segment> read_axis(_parser, _parser.read_int(_par[1]));
    read_axis.read_line_segment();
    Read_entity<Entity<1>> read_generatrix(_parser, _parser.read_int(_par[2]));
    read_generatrix.read_curve();
    Entity<1>* generatrix = read_generatrix.get().release();
    generatrix->scale = 1.;
    Line_segment axis = *read_axis.get();
    axis.scale = 1.;
    auto surf = new Revolution_surface {
      generatrix,
      axis,
      _parser.read_float(_par[3]),
      _parser.read_float(_par[4]),
    };
    _ptr.reset(surf);
  }

  void read_surface() {
    read_plane();
    read_revolution_surface();
    HEXED_ASSERT(_ptr, "failed to read surface from entity #" + std::to_string(_ent_num));
  }

  void read_disc_curve() {
    Read_entity<Entity<1>> curve(_parser, _dir);
    curve.read_curve();
    auto ptr = curve.get();
    _ptr.reset(new Composite_curve);
    if (ptr) {
      _ptr->emplace_back(ptr.release());
      _ptr->back()->scale = 1.;
    } else if (curve.entity_number() == 116 || curve.entity_number() == 132) {
    } else if (curve.entity_number() == 102) {
      for (int i_curve = 0; i_curve < _parser.read_int(_par[1]); ++i_curve) {
        Read_entity<T> sub_curve(_parser, _parser.read_int(_par[2 + i_curve]));
        sub_curve.read_disc_curve();
        for (auto& c : *sub_curve._ptr) _ptr->emplace_back(c.release());
      }
    } else HEXED_THROW("failed to read curve from entity #" + std::to_string(curve.entity_number()));
  }

  void read_trimmed_surface() {
    if (_ptr || _ent_num != 144) return;
    Read_entity<Entity<2>> read_surf(_parser, _parser.read_int(_par[1]));
    read_surf.read_surface();
    if (read_surf.entity_number() == 108) return;
    _ptr.reset(new Trimmed_surface());
    _ptr->surface.reset(read_surf.get().release());
    _ptr->surface->scale = 1.;
    HEXED_ASSERT(_parser.read_int(_par[2]) == 1, "using outer boundary as boundary curve is not implemented",
                 assert::Not_implemented_error);
    double sz = 1./_ptr->n_div;
    std::vector<std::vector<Mat<2>>> curves;
    Mat<2, 2> bounds;
    bounds << huge, -huge, huge, -huge;
    for (Int i_curve = 0; i_curve < 1 + _parser.read_int(_par[3]); ++i_curve) {
      Read_entity<T> on_surf(_parser, _parser.read_int(_par[4 + i_curve]));
      HEXED_ASSERT(on_surf.entity_number() == 142, "boundary must be a curve on a surface");
      int model_curve = _parser.read_int(on_surf._par[4]);
      std::vector<Mat<2>> nodes;
      if ((_parser.read_int(on_surf._par[5]) == 1 || model_curve == 0) && _parser.read_int(on_surf._par[3]) != 0) {
        HEXED_THROW("parameter-space curves are not implemented", assert::Not_implemented_error);
      } else {
        HEXED_ASSERT(model_curve != 0,
                     "at least one of parameter-space and model-space curve pointers must be defined");
        Read_entity<Composite_curve> curve(_parser, model_curve);
        curve.read_disc_curve();
        _ptr->curves.emplace_back(curve.get().release());
        for (auto& c : *_ptr->curves.back()) {
          for (Int i_div = 0; i_div < _ptr->n_div + 1; ++i_div) {
            Mat<3> pt = c->point(Mat<1>{i_div*sz});
            Mat<2> params = _ptr->surface->nearest_params(pt, [](Mat<2>){return true;}).params;
            nodes.push_back(params);
          }
        }
      }
      if (!nodes.empty()) {
        for (Mat<2> node : nodes) {
          bounds(all, 0) = bounds(all, 0).cwiseMin(node);
          bounds(all, 1) = bounds(all, 1).cwiseMax(node);
        }
      }
      curves.push_back(std::move(nodes));
    }
    bounds(all, 1) = bounds(all, 1).cwiseMax(bounds(all, 0) + Mat<2>{sz, sz});
    _ptr->surface->reparameterize(bounds);
    for (auto& nodes : curves) if (!nodes.empty()) {
      for (int i_node = 0; i_node < Int(nodes.size()); ++i_node) {
        nodes[i_node] = (nodes[i_node] - bounds(all, 0)).cwiseQuotient(bounds(all, 1) - bounds(all, 0));
      }
      Int n_nodes = nodes.size();
      for (Int i_node = n_nodes, changed = false; (i_node < 2*n_nodes) || changed; ++i_node, changed = false) {
        for (int i_dim = 0; i_dim < 2; ++i_dim) {
          for (int sign : {-1, 1}) {
            if (std::abs(nodes[i_node%n_nodes](i_dim) + sign - nodes[(i_node - 1)%n_nodes](i_dim)) <
                std::abs(nodes[i_node%n_nodes](i_dim)        - nodes[(i_node - 1)%n_nodes](i_dim))) {
              nodes[i_node%nodes.size()](i_dim) += sign;
              changed = true;
            }
          }
        }
      }
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        Int n_less = 0;
        Int n_greater = 0;
        for (Mat<2> node : nodes) {
          n_less += node(i_dim) < -sz;
          n_greater += node(i_dim) > 1 + sz;
        }
        int sign = (n_less > n_nodes/2) - (n_greater > n_nodes/2);
        for (Mat<2>& node : nodes) node(i_dim) += sign;
      }
      std::vector<Int> abscissa;
      std::vector<double> ordinate;
      Mat<2> prev_params {-1., 0.};
      Int prev_absc = -1;
      if (n_nodes) for (Int i_node = 0; i_node < n_nodes + 1; ++i_node) {
        Mat<2> params;
        if (i_node == n_nodes && !abscissa.empty()) params << abscissa.front(), ordinate.front();
        else {
          params = nodes[i_node];
          params(0) = std::max(0., std::min(1., params(0)))*_ptr->n_div;
        }
        if (prev_params(0) < -.1) prev_params = params;
        while (floor(params(0)) > floor(prev_params(0)) || ceil(params(0)) < ceil(prev_params(0))) {
          prev_absc = floor(params(0)) > floor(prev_params(0)) ? floor(prev_params(0)) + 1
                                                               : ceil (prev_params(0)) - 1;
          double denom = params(0) - prev_params(0);
          if (std::abs(denom) < sz) prev_params(1) = params(1);
          else prev_params(1) += (prev_absc - prev_params(0))*(params(1) - prev_params(1))/denom;
          prev_params(0) = prev_absc;
          abscissa.push_back(prev_absc);
          ordinate.push_back(prev_params(1));
        }
      }
      if (!abscissa.empty()) {
        HEXED_ASSERT(abscissa.front() == abscissa.back(), "parametric representation is not closed");
        for (Int i_segment = 0; i_segment < Int(abscissa.size()) - 1; ++i_segment) {
          HEXED_ASSERT(std::abs(abscissa[i_segment] - abscissa[i_segment + 1]) == 1, "invalid step size");
          bool reverse = abscissa[i_segment] > abscissa[i_segment + 1];
          HEXED_ASSERT(0 <= abscissa[i_segment + reverse] && abscissa[i_segment + reverse] < _ptr->n_div,
                       format_str(200, "segment index %i out of bounds", abscissa[i_segment + reverse]));
          _ptr->parametric_segments[abscissa[i_segment + reverse]].push_back({ordinate[i_segment +  reverse],
                                                                              ordinate[i_segment + !reverse]});
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
      _ptr->scale = _unit(units[unit_flag - 1]);
    }
    return std::unique_ptr<T>{_ptr.release()};
  }

  Int entity_number() const {return _ent_num;}

  private:
  double _unit(double unit) const {return unit;}
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
    if (ts) {
      for (auto& composite : ts->curves) {
        for (auto& curve : *composite) {
          Array<double> nodes({ts->n_div + 1, 3});
          for (Int i_node = 0; i_node < ts->n_div + 1; ++i_node) {
            nodes(i_node).vector() = ts->_convert(curve->point(Mat<1>{i_node/double(ts->n_div)}));
          }
          _edges.emplace_back(nodes);
        }
      }
      _surfaces.emplace_back(ts.release());
    }
  }
}

Nearest_point<dyn> Geom::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  Nearest_point<dyn> nearest(point, distance_guess);
  for (auto& surf : _surfaces) nearest.merge(Mat<>{surf->nearest_point(point, distance_guess)});
  if ((!nearest.empty() && std::sqrt(nearest.dist_squared()) < distance_guess) || distance_guess >= max_distance) {
    std::cout << "| " << nearest.dist_squared() << std::endl;
    return nearest;
  }
  std::cout << distance_guess << " ";
  return nearest_point(point, max_distance, distance_guess*2);
}

void Geom::visualize(std::string file_name) const {
  int n = 101;
  for (Int i_edge = 0; i_edge < Int(_edges.size()); ++i_edge) {
    _edges[i_edge].visualize("default", file_name + "edge" + std::to_string(i_edge));
  }
  {
    auto vis = Visualizer::create("default", 3, 1, "parametric_edges", {}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      for (Int i_div = 0; i_div < s->n_div; ++i_div) {
        for (auto& seg : s->parametric_segments[i_div]) {
          Array<double> coords({3, 2});
          for (int i_node = 0; i_node < 2; ++i_node) {
            Mat<2> params;
            params(0) = (i_div + i_node)/double(s->n_div);
            params(1) = seg(i_node);
            Mat<3> point = s->point(params);
            for (int i_dim = 0; i_dim < 3; ++i_dim) coords(i_dim)[i_node] = point(i_dim);
          }
          vis->write_block(coords, Array<double>({0, 2}));
        }
      }
    }
  }
  {
    auto vis = Visualizer::create("default", 3, 2, "surfaces", {"inside"}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      Array<double> discrete({3, n, n});
      Array<double> inside({1, n, n});
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          Mat<2> params {i/(n - 1.), j/(n - 1.)};
          Mat<3> p = s->point(params);
          for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)(i)[j] = p(i_dim);
          inside(0)(i)[j] = s->inside(params);
        }
      }
      vis->write_block(discrete, inside);
    }
  }
  {
    auto vis = Visualizer::create("default", 3, 3, "distance", {"distance"}, 0., Visualizer::block);
    Array<double> coords({3, n, n, n});
    Array<double> dist({1, n, n, n});
    #pragma omp parallel for
    for (int i = 0; i < math::pow(n, 3); ++i) {
      Mat<3> p;
      for (int i_dim = 0; i_dim < 3; ++i_dim) p(i_dim) = coords(i_dim)[i] = 1./n*(i/math::pow(n, 2 - i_dim)%n) - .1;
      double min_dist = huge;
      for (auto& s : _surfaces) min_dist = std::min(min_dist, (p - s->nearest_point(p, huge)).norm());
      dist[i] = min_dist;
    }
    vis->write_block(coords, dist);
  }
}

}
