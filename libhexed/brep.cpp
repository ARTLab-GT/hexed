#include <optional>
#include <hexed/brep.hpp>
#include <hexed/Iges_parser.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/utils.hpp>
#include <hexed/constants.hpp>

namespace hexed::brep {

//! \cond

Line_segment::Line_segment(Mat<3, 2> endpoints)
: _endpoints{endpoints}
, _length{(_endpoints(all, 1) - _endpoints(all, 0)).norm()}
{}

Parametric<1>::Nearest_parameters Line_segment::nearest_params(Mat<3> p, Constraint is_feasible,
                                                               double max_distance) const {
  Mat<3> diff = _endpoints(all, 1) - _endpoints(all, 0);
  Mat<1> params {std::max(0., std::min(1., (p - _endpoints(all, 0)).dot(diff)/diff.squaredNorm()))};
  return {params, is_feasible(params)};
}

std::vector<Parametric<1>::Intersection_parameters> Line_segment::intersection_params(Mat<3, 2> points) const {
  Mat<3, 2> lhs;
  lhs(all, 0) = _endpoints(all, 1) - _endpoints(all, 0);
  lhs(all, 1) = points(all, 0) - points(all, 1);
  lhs(2, all).setZero();
  Mat<3> rhs = points(all, 0) - _endpoints(all, 0);
  rhs(2) = 0;
  Mat<2> soln = lhs.householderQr().solve(rhs);
  if (0 <= soln(0) && soln(0) <= 1) return {{Mat<1>{soln(0)}, soln(1)}};
  return {};
}

Circular_arc::Circular_arc(Mat<3> center, double radius, double start_angle, double end_angle)
: _center{center}
, _radius{radius}
, _start_angle{start_angle}
, _end_angle{end_angle}
{
  HEXED_ASSERT(_end_angle - _start_angle > 0, "end angle must be greater than start angle");
}

double limited_angle(double angle, double start, double end) {
  if (math::angle_diff(angle, start) > end - start) {
    return (math::angle_diff(angle, end) > math::angle_diff(start, angle)) ? start : end;
  }
  return angle;
}

Mat<3> Circular_arc::point(Mat<1> params) const {
  double angle = _start_angle + params(0)*(_end_angle - _start_angle);
  return _center + _radius*Mat<3>{std::cos(angle), std::sin(angle), 0.};
}

Parametric<1>::Nearest_parameters Circular_arc::nearest_params(Mat<3> p, Constraint is_feasible,
                                                               double max_distance) const {
  double angle = limited_angle(std::atan2(p(1) - _center(1), p(0) - _center(0)), _start_angle, _end_angle);
  Mat<1> params {math::angle_diff(angle, _start_angle)/(_end_angle - _start_angle)};
  return {params, is_feasible(params)};
}

std::vector<Parametric<1>::Intersection_parameters> Circular_arc::intersection_params(Mat<3, 2> points) const {
  Mat<3> start = points(all, 0) - _center;
  Mat<3> diff = points(all, 1) - points(all, 0);
  // radius^2 = coefs[0] + coefs[1]*interp_coef + coefs[2]*interp_coef^2
  double coefs [3] {
    start(0)*start(0) + start(1)*start(1),
    2*(start(0)*diff(0) + start(1)*diff(1)),
    diff(0)*diff(0) + diff(1)*diff(1),
  };
  double descrim = coefs[1]*coefs[1] - 4*coefs[2]*(coefs[0] - _radius*_radius);
  std::vector<Parametric<1>::Intersection_parameters> sects;
  // skip the single intersection case---no one's actually going to care
  if (descrim > 0) {
    double root = std::sqrt(descrim);
    for (int sign : {-1, 1}) {
      double interp = (-coefs[1] + sign*root)/(2*coefs[2]);
      Mat<3> sect_point = start + interp*diff;
      double param = math::angle_diff(std::atan2(sect_point(1), sect_point(0)), _start_angle)
                     /(_end_angle - _start_angle);
      if (0 <= param && param <= 1) sects.push_back({Mat<1>{param}, interp});
    }
  }
  return sects;
}

Parametric<2>::Nearest_parameters Plane::nearest_params(Mat<3> p, Constraint is_feasible,
                                                        double max_distance) const {
  Mat<2> params = _vecs.colPivHouseholderQr().solve(p - _origin);
  return {params, is_feasible(params)};
}

std::vector<Parametric<2>::Intersection_parameters> Plane::intersection_params(Mat<3, 2> endpoints) const {
  Mat<3, 3> lhs;
  lhs(all, Eigen::seqN(0, 2)) = _vecs;
  lhs(all, 2) = endpoints(all, 0) - endpoints(all, 1);
  Mat<3> rhs = endpoints(all, 0) - _origin;
  auto fact = lhs.fullPivHouseholderQr();
  if (!fact.isInvertible()) return {};
  Mat<3> soln = fact.solve(rhs);
  for (int i = 0; i < 2; ++i) if (!(soln(i) >= 0 && soln(i) <= 1)) return {};
  return {{{soln(0), soln(1)}, soln(2)}};
}

Mat<2, 2> Plane::reparameterize(Mat<2, 2> bounds) {
  _origin = point(bounds(all, 0));
  _vecs = _vecs*(bounds(all, 1) - bounds(all, 0)).asDiagonal();
  return bounds;
}

// discretizes a curve into `n_div` polygonal segments and returns their `n_div + 1` endpoints
Array<double> discretize(Parametric<1>& curve, Int n_div) {
  Array<double> nodes({n_div + 1, 3});
  for (Int i_node = 0; i_node < n_div + 1; ++i_node) {
    nodes(i_node).vector() = curve.point(Mat<1>{i_node/double(n_div)});
  }
  return nodes;
}

Revolution_surface::Revolution_surface(Parametric<1>* g, Line_segment ax, Int n_div, double sa, double ea)
: _generatrix{g}
, _axis{ax}
, _unit_axis{(_axis.point(Mat<1>{1.}) - _axis.point(Mat<1>{0.})).normalized()}
, _n_div{n_div}
, _start_angle{sa}
, _end_angle{ea}
, _tree(discretize(*_generatrix, _n_div), 4)
{
  HEXED_ASSERT(_end_angle - _start_angle > 0, "end angle must be greater than start angle");
}

Mat<3> Revolution_surface::rotate(Mat<3> p, double angle) const {
  // compute the displacement relative to the first endpoint
  Mat<3> origin = _axis.point(Mat<1>{0.});
  p -= origin;
  // separate into components parallel and orthogonal to the axis
  Mat<3> ax_vec = (_axis.point(Mat<1>{1.}) - origin).normalized();
  Mat<3> axial_component = p.dot(ax_vec)*ax_vec;
  Mat<3> radial_component = p - axial_component;
  // rotate the radial (orthogonal to axis) component and recombine
  return std::cos(angle)*radial_component + std::sin(angle)*ax_vec.cross(radial_component)
         + axial_component + origin;
};

Mat<3> Revolution_surface::point(Mat<2> params) const {
  Mat<3> p = _generatrix->point(Mat<1>{params(0)});
  double angle = _start_angle + params(1)*(_end_angle - _start_angle);
  return rotate(p, angle);
}

// given a point `arc_point` which is nominally on the genratrix, compute the rotation angle of the nearest point
// on the arc of points on the surface obtained by rotating `arc_point`
double Revolution_surface::_unlimited_best_angle(Mat<3> arc_point, Mat<3> radius) const {
  arc_point -= _axis.point(Mat<1>{0.});
  Mat<3> arc_radius = (arc_point - arc_point.dot(_unit_axis)*_unit_axis).normalized();
  return std::atan2(arc_radius.cross(radius).dot(_unit_axis), arc_radius.dot(radius));
}

double Revolution_surface::_limited_best_angle(Mat<3> arc_point, Mat<3> radius) const {
  return limited_angle(_unlimited_best_angle(arc_point, radius), _start_angle, _end_angle);
}

// given a point `arc_point` which is nominally on the genratrix, compute the nearst point
// on the arc of points on the surface obtained by rotating `arc_point`
Mat<3> Revolution_surface::_best_point(Mat<3> arc_point, Mat<3> radius) const {
  double angle = _limited_best_angle(arc_point, radius);
  return rotate(arc_point, angle);
}

// helper class that does the real work of nearest point calculations
class Revolution_surface::_Find_nearest {
  public:
  struct Candidate {
    Nearest_parameters np;
    double dist;
  };
  const Revolution_surface& surf;
  Constraint is_feasible;
  Mat<3> point; // point we want to compute the nearest point to
  Mat<3> from_start; // displacement of `point` relative to first axis endpoint
  Mat<3> radius; // component of `from_start` orthogonal to axis
  Candidate cand;
  _Find_nearest(const Revolution_surface& s, Mat<3> p, Parametric<2>::Constraint is_f, double max_distance)
  : surf{s}
  , is_feasible{is_f}
  , point{p}
  , from_start{p - surf._axis.point(Mat<1>{0.})}
  , radius{(from_start - from_start.dot(surf._unit_axis)*surf._unit_axis).normalized()}
  // initialize `cand` to a non-point but set the distance to `max_distance`
  // so that any farther candidates will be ignored
  , cand{{Mat<2>{std::nan(""), std::nan("")}, false}, max_distance}
  {}
  // return whichever of `c0` and `c1` is a better candidate
  Candidate merge(Candidate c0, Candidate c1) {
    if (c1.dist < c0.dist && c1.np.is_feasible) return c1;
    return c0;
  }
  // search the part of the surface subtended by the segment of the generatrix covered by the nodes of `segment`
  // for the nearest point.
  // Any candidates for the nearest point in this area will be `merge`d with `cand`.
  void find(const Tree_curve::Segment& segment) {
    if (segment.segments.size()) {
      // if this segment is not a leaf, recursively search its child segments
      double dist [2];
      for (int i_segment = 0; i_segment < 2; ++i_segment) {
        dist[i_segment] = (surf._best_point(segment.segments[i_segment].center, radius) - point).norm();
      }
      // do the closer segment first in hopes that we can find a point close enough
      // to justify skipping the farther one
      for (bool i_segment : {dist[1] < dist[0], !(dist[1] < dist[0])}) {
        if (dist[i_segment] - segment.segments[i_segment].radius < cand.dist) find(segment.segments[i_segment]);
      }
    } else {
      // this segment is a leaf, so search its nodes
      Int n_nodes = segment.nodes.shape()[0];
      for (Int i_node = 0; i_node < n_nodes; ++i_node) {
        // find the parameters of the neares point on the arc subtended by this node
        Candidate c;
        Mat<3> node = segment.nodes(i_node).vector();
        c.np.params(0) = double(segment.nodes_start + i_node)/surf._n_div;
        double angle = surf._limited_best_angle(node, radius);
        c.np.params(1) = math::angle_diff(angle, surf._start_angle)/(surf._end_angle - surf._start_angle);
        // check if the computed nearest parameters are feasible
        c.np.is_feasible = is_feasible(c.np.params);
        c.dist = (surf.rotate(node, angle) - point).norm();
        // merge candidates
        cand = merge(cand, c);
      }
    }
  }
};

class Revolution_surface::_Find_intersects {
  public:
  const Revolution_surface& surf;
  std::vector<Intersection_parameters> intersects;
  struct Coefs {
    Mat<3> radial;
    Mat<3> axial;
    Coefs(const Revolution_surface& surf, Mat<3, 2> points) {
      Mat<3> r_start = points(all, 0) - surf._axis.point(Mat<1>{0.});
      axial[0] = r_start.dot(surf._unit_axis);
      r_start -= axial[0]*surf._unit_axis;
      Mat<3> r_diff = points(all, 1) - points(all, 0);
      axial[1] = r_diff.dot(surf._unit_axis);
      r_diff -= axial[1]*surf._unit_axis;
      radial[0] = r_start.squaredNorm();
      radial[1] = 2*r_start.dot(r_diff);
      radial[2] = r_diff.squaredNorm();
    }
    Coefs() = default;
  };
  Mat<3, 2> points;
  Coefs points_coefs;
  Int total_nodes;
  _Find_intersects(const Revolution_surface& s, Mat<3, 2> p)
  : surf{s}, points{p}, points_coefs{s, points}, total_nodes(surf._tree.n_points())
  {}
  void find(const Tree_curve::Segment& segment) {
    if (segment.segments.size()) {
      for (const Tree_curve::Segment& segment : segment.segments) {
        Mat<3> center = segment.center - surf._axis.point(Mat<1>{0.});
        double axial = center.dot(surf._unit_axis);
        double radius = (center - axial*surf._unit_axis).norm();
        bool could_intersect;
        if (points_coefs.radial[2] < math::pow(points_coefs.axial[1], 2)) {
          double center_param = (axial - points_coefs.axial[0])/points_coefs.axial[1];
          double param_diff = segment.radius/points_coefs.axial[1];
          double rsq = points_coefs.radial[0]
                       + (points_coefs.radial[1] + points_coefs.radial[2]*center_param)*center_param;
          double rsq_deriv = points_coefs.radial[1] + 2*points_coefs.radial[2]*center_param;
          double rsq_uncert = (std::abs(rsq_deriv) + points_coefs.radial[2]*param_diff)*param_diff;
          could_intersect = std::abs(radius*radius - rsq) <
                            segment.radius*segment.radius + 2*std::abs(radius*segment.radius) + rsq_uncert;
        } else {
          could_intersect = true;
        }
        if (could_intersect) find(segment);
      }
    } else {
      Int n_nodes = segment.nodes.shape()[0];
      Coefs coefs [2];
      coefs[0] = points_coefs;
      for (Int i_node = 0; i_node < n_nodes - 1; ++i_node) {
        Mat<3, 2> gener_points;
        for (int col = 0; col < 2; ++col) gener_points(all, col) = segment.nodes(i_node + col).vector();
        coefs[1] = Coefs(surf, gener_points);
        int i_transform = std::abs(coefs[0].axial[1]) < std::abs(coefs[1].axial[1]);
        Mat<2> transform {
          (coefs[!i_transform].axial[0] - coefs[i_transform].axial[0])/coefs[i_transform].axial[1],
          coefs[!i_transform].axial[1]/coefs[i_transform].axial[1],
        };
        Mat<3> quad_coefs = coefs[i_transform].radial;
        quad_coefs(0) += transform(0)*(quad_coefs(1) + transform(0)*quad_coefs(2));
        quad_coefs(1) += 2*transform(0)*quad_coefs(2);
        for (int pow = 0; pow < 3; ++pow) quad_coefs(pow) *= math::pow(transform(1), pow);
        quad_coefs -= coefs[!i_transform].radial;
        if (quad_coefs(2) != 0 && std::isfinite(transform(1))) {
          double descrim = quad_coefs(1)*quad_coefs(1) - 4*quad_coefs(2)*quad_coefs(0);
          if (descrim > 0) {
            for (int sign : {-1, 1}) {
              double soln = (-quad_coefs(1) + sign*std::sqrt(descrim))/(2*quad_coefs(2));
              double node_interp = i_transform ? transform(0) + transform(1)*soln : soln;
              if (-1e-3 <= node_interp && node_interp <= 1 + 1e-3) {
                double interp_coef = i_transform ? soln : transform(0) + transform(1)*soln;
                double param0 = (segment.nodes_start + i_node + node_interp)/(total_nodes - 1);
                Mat<3> gen_point = gener_points*Mat<2>{1. - node_interp, node_interp};
                Mat<3> radius = points*Mat<2>{1. - interp_coef, interp_coef} - surf._axis.point(Mat<1>{0.});
                radius -= radius.dot(surf._unit_axis)*surf._unit_axis;
                double angle = surf._unlimited_best_angle(gen_point, radius.normalized());
                double param1 = math::angle_diff(angle, surf._start_angle)/(surf._end_angle - surf._start_angle);
                if (param1 <= 1) intersects.push_back({{param0, param1}, interp_coef});
              }
            }
          }
        }
      }
    }
  }
};

Parametric<2>::Nearest_parameters Revolution_surface::nearest_params(Mat<3> p, Constraint is_feasible,
                                                                     double max_distance) const {
  _Find_nearest finder(*this, p, is_feasible, max_distance);
  finder.find(_tree.root());
  return finder.cand.np;
}

std::vector<Parametric<2>::Intersection_parameters> Revolution_surface::intersection_params(Mat<3, 2> pnts) const {
  _Find_intersects finder(*this, pnts);
  finder.find(_tree.root());
  return finder.intersects;
}

Coordinate_change::Coordinate_change(Mat<3> translate, Mat<3, 3> transform)
: _translate{translate}, _transform{transform}, _inv{transform.inverse()}
{}

Coordinate_change Coordinate_change::operator()(Coordinate_change that) const {
  return {_translate + _transform*that._translate, _transform*that._transform};
}

Trimmed_surface::Trimmed_surface(Parametric<2>* surface, std::vector<Composite_curve>&& curves, Int n_div)
: _n_div{n_div}, _sz{1./_n_div}, _surf{surface}
{
  // discretize curves into polygonal segments in parameter space
  std::vector<std::vector<Mat<2>>> discrete_curves;
  for (auto& composite : curves) {
    discrete_curves.emplace_back();
    auto& param_nodes = discrete_curves.back();
    for (auto& curve : composite) {
      Array<double> phys_nodes({_n_div + 1, 3});
      Mat<3> start = curve->point(Mat<1>{0.});
      double mean_squared_dist = 0;
      for (Int i_node = 0; i_node < _n_div + 1; ++i_node) {
        Mat<3> pt = curve->point(Mat<1>{i_node*_sz});
        mean_squared_dist += (pt - start).squaredNorm();
        phys_nodes(i_node).vector() = pt;
        Mat<2> params = _surf->nearest_params(pt, [](Mat<2>){return true;}, default_max_dist).params;
        param_nodes.push_back(params);
      }
      mean_squared_dist /= n_div + 1;
      if ((curve->point(Mat<1>{1.}) - start).squaredNorm() < .1*mean_squared_dist) {
        _curves.emplace_back(phys_nodes(0, n_div/2 + 1).copy(), 4);
        _curves.emplace_back(phys_nodes(n_div/2, n_div + 1).copy(), 4);
      } else {
        _curves.emplace_back(phys_nodes.copy(), 4);
      }
    }
    if (!param_nodes.empty()) param_nodes.push_back(param_nodes.front());
  }
  // initialize parameter-space curves with discretiation
  initialize(discrete_curves);
}

void Trimmed_surface::initialize(std::vector<std::vector<Mat<2>>>& curves) {
  // compute bounds of discrete nodes in parameter space
  Mat<2, 2> bounds;
  bounds << huge, -huge, huge, -huge;
  bool set = false;
  for (auto& curve : curves) {
    set = set || !curve.empty();
    for (Mat<2> node : curve) {
      bounds(all, 0) = bounds(all, 0).cwiseMin(node);
      bounds(all, 1) = bounds(all, 1).cwiseMax(node);
    }
  }
  if (set) bounds(all, 1) = bounds(all, 1).cwiseMax(bounds(all, 0) + Mat<2>{_sz, _sz});
  else bounds << 0, 1, 0, 1;
  // reparameterize surface to contain bounds
  bounds = _surf->reparameterize(bounds);
  _param_segments.resize(_n_div);
  for (auto& nodes : curves) if (!nodes.empty()) {
    Int n_nodes = nodes.size();
    // apply reparameterization to nodes
    for (int i_node = 0; i_node < Int(nodes.size()); ++i_node) {
      nodes[i_node] = (nodes[i_node] - bounds(all, 0)).cwiseQuotient(bounds(all, 1) - bounds(all, 0));
    }
    // correct periodic seam errors
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
    // the seam correction process can end up shifting the entire curve loop to [-1, 0] or [1, 2],
    // so this loop shifts it again to keep the largest possible number of nodes in [0, 1]
    for (int i_dim = 0; i_dim < 2; ++i_dim) {
      Int n_less = 0;
      Int n_greater = 0;
      for (Mat<2> node : nodes) {
        n_less += node(i_dim) < -_sz;
        n_greater += node(i_dim) > 1 + _sz;
      }
      int sign = (n_less > n_nodes/2) - (n_greater > n_nodes/2);
      for (Mat<2>& node : nodes) node(i_dim) += sign;
    }
    // compute parametric segments
    std::vector<Int> abscissa;
    std::vector<double> ordinate;
    Mat<2> prev_params {-1., 0.};
    Int prev_absc = -1;
    if (n_nodes) for (Int i_node = 0; i_node <= n_nodes; ++i_node) {
      Mat<2> params;
      if (i_node == n_nodes && !abscissa.empty()) params << abscissa.front(), ordinate.front();
      else {
        params = nodes[i_node];
        params(0) = std::max(0., std::min(1., params(0)))*_n_div;
      }
      if (prev_params(0) < -.1) prev_params = params;
      while (floor(params(0)) > floor(prev_params(0)) || ceil(params(0)) < ceil(prev_params(0))) {
        prev_absc = floor(params(0)) > floor(prev_params(0)) ? floor(prev_params(0)) + 1
                                                             : ceil (prev_params(0)) - 1;
        double denom = params(0) - prev_params(0);
        if (std::abs(denom) < _sz) prev_params(1) = params(1);
        else prev_params(1) += (prev_absc - prev_params(0))*(params(1) - prev_params(1))/denom;
        prev_params(0) = prev_absc;
        abscissa.push_back(prev_absc);
        ordinate.push_back(prev_params(1));
      }
    }
    // sort segments into bins of specified `param(0)`
    if (!abscissa.empty()) {
      HEXED_ASSERT(abscissa.front() == abscissa.back(), "parametric representation is not closed");
      for (Int i_segment = 0; i_segment < Int(abscissa.size()) - 1; ++i_segment) {
        HEXED_ASSERT(std::abs(abscissa[i_segment] - abscissa[i_segment + 1]) == 1, "invalid step size");
        bool reverse = abscissa[i_segment] > abscissa[i_segment + 1];
        HEXED_ASSERT(0 <= abscissa[i_segment + reverse] && abscissa[i_segment + reverse] < _n_div,
                     format_str(200, "segment index %i out of bounds", abscissa[i_segment + reverse]));
        _param_segments[abscissa[i_segment + reverse]].push_back({ordinate[i_segment +  reverse],
                                                                  ordinate[i_segment + !reverse]});
      }
    }
  }
}

bool Trimmed_surface::is_inside(Mat<2> params) const {
  // count the number of segmements intersected by a ray in the positive `params(1)` direction
  Int i_seg = floor(params(0)*_n_div);
  if (i_seg < 0 || i_seg > _n_div) return false;
  if (i_seg == _n_div) i_seg = _n_div - 1;
  Int n_intersections = 0;
  for (Mat<2> seg : _param_segments[i_seg]) n_intersections += params(1) < seg(0) + (params(0)*_n_div - i_seg)*(seg(1) - seg(0));
  // the point is inside iff the number of intesections is odd
  return n_intersections%2;
}

next::Sequence<const Tree_curve&> Trimmed_surface::curves() const {
  return next::Sequence<const Tree_curve&>::vector_view(_curves);
}

Nearest_point<3> Trimmed_surface::nearest_point(Mat<3> point, double max_dist) const {
  Nearest_point<3> nearest(point, max_dist);
  // first find the nearest point in the surface interior, if any
  auto params = _surf->nearest_params(point, [this](Mat<2> params){return is_inside(params);}, max_dist);
  if (params.is_feasible) nearest.merge(_surf->point(params.params));
  // then check the nearest point on all the boundary curves
  if (nearest.empty() || _surf->must_check_boundary()) {
    for (auto& curve : _curves) {
      auto index = curve.nearest_point(point, 1.01*std::sqrt(nearest.dist_squared()));
      if (index.index > -1) nearest.merge(curve.nodes()(index.index).vector());
    }
  }
  return nearest;
}

std::vector<double> Trimmed_surface::intersections(Mat<3, 2> endpoints) const {
  std::vector<double> sects;
  auto sect_params = _surf->intersection_params(endpoints);
  for (auto params : sect_params) {
    if (is_inside(params.params)) sects.push_back(params.interp_coef);
  }
  return sects;
}

// helper class to read an entity from an IGES file
class Read_entity {
  public:
  // Initializes the to read the entity whose directory entry starts on line `line`.
  // If `change_to` is provided, it will be postcomposed with any coordinate transformation specified by the entity.
  Read_entity(const Iges_parser& parser, Int line, Int n_div, Coordinate_change change_to = {})
  : Read_entity{parser, parser.entry(Iges_parser::directory, line), n_div, change_to}
  {}
  // initializes the reader to read the entity whose directory entry is given literally by `dir`
  Read_entity(const Iges_parser& parser, const std::vector<std::string>& dir, Int n_div,
              Coordinate_change change_to = {})
  : _parser{parser}
  , _dir{dir}
  , _par{_parser.entry(Iges_parser::parameter, _parser.read_int(_dir[1]))}
  , _ent_num{_parser.read_int(_dir[0])}
  , _n_div{n_div}
  {
    // compute the coordinate transformation of this entity
    Int line = _parser.read_int(_dir[6]);
    if (!line) _coords = change_to; // IGES entities can have a null pointer for indicating the identity transformation
    else {
      Read_entity reader(_parser, line, _n_div);
      _coords = change_to(reader._read_coord());
    }
    // determine the units of the file
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
    _unit = units[unit_flag - 1];
  }

  template <typename T, typename U>
  static void merge(std::unique_ptr<T>& ptr0, std::unique_ptr<U>&& ptr1) {
    if (!ptr0 && ptr1) ptr0.reset(ptr1.release());
  }

  /* The following `read_` functions all attempt to read a specific entity type.
   * If the entry `this` is reading is indeed that entity type, the corresponding `hexed::brep` object
   * will be constructed and a `unique_ptr` to it will be returned.
   * Otherwise, an empty `unique_ptr` will be returned.
   * Thus, if you don't know what type of entity you're trying to read, you can just spam all the `read_` functions
   * and see which one gives you a non-empty value.
   */

  std::unique_ptr<Line_segment> read_line_segment() {
    if (_ent_num != 110) return {};
    Mat<3, 2> endpts;
    for (int i = 0; i < 6; ++i) endpts(i) = _unit*_parser.read_float(_par[1 + i]);
    for (int i = 0; i < 2; ++i) endpts(all, i) = _coords.to_model(endpts(all, i));
    return std::make_unique<Line_segment>(endpts);
  }

  std::unique_ptr<Parametric<1>> read_circular_arc() {
    if (_ent_num != 100) return {};
    std::vector<double> values;
    for (std::string s : _par) values.push_back(_parser.read_float(s));
    Mat<3> center {values[2], values[3], values[1]};
    center *= _unit;
    double radius = _unit*std::sqrt(.5*(values[4]*values[4] + values[5]*values[5]
                                    + values[6]*values[6] + values[7]*values[7]));
    double start_angle = std::atan2(values[5], values[4]);
    double end_angle = std::atan2(values[7], values[6]);
    while (end_angle < start_angle + 1e-10) end_angle += 2*constants::pi;
    std::unique_ptr<Parametric<1>> arc (new Circular_arc {
      center, radius,
      start_angle,
      end_angle,
    });
    return std::unique_ptr<Parametric<1>>(new Transformed<1>(arc.release(), _coords));
  }

  std::unique_ptr<Plane> read_plane() const {
    if (_ent_num != 108) return {};
    HEXED_ASSERT(_parser.read_int(_par[5]) == 0, "bounded planes are not implemented", assert::Not_implemented_error);
    double coefs [4];
    for (int i_coef = 0; i_coef < 4; ++i_coef) coefs[i_coef] = _parser.read_float(_par[i_coef + 1]);
    coefs[3] *= _unit;
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
    return std::make_unique<Plane>(_coords.to_model(origin), _coords.transform()*vecs);
  }

  std::unique_ptr<Revolution_surface> read_revolution_surface() const {
    if (_ent_num != 120) return {};
    Read_entity read_axis(_parser, _parser.read_int(_par[1]), _n_div, _coords);;
    auto axis = read_axis.read_line_segment();
    Read_entity read_generatrix(_parser, _parser.read_int(_par[2]), _n_div, _coords);
    auto generatrix = read_generatrix.read_curve();
    return std::make_unique<Revolution_surface>(generatrix.release(), *axis, _n_div,
                                                _parser.read_float(_par[3]), _parser.read_float(_par[4]));
  }

  // Attempts to read any of the entities that derive from `Parametric<1>`.
  // Iff `required == true`, throws on failure.
  std::unique_ptr<Parametric<1>> read_curve(bool required = true) {
    std::unique_ptr<Parametric<1>> ptr;
    merge(ptr, read_line_segment());
    merge(ptr, read_circular_arc());
    HEXED_ASSERT(
      !required || ptr,
      "Curve entity #" + std::to_string(_ent_num) + " is not implemented.",
      assert::Not_implemented_error
    );
    return ptr;
  }

  // attempts to read a composite curve (collection of curves that share endpoints) and throws on failure
  Composite_curve read_composite_curve() {
    auto ptr = read_curve(false);
    Composite_curve comp;
    if (ptr) {
      comp.emplace_back(ptr.release());
    } else if (_ent_num == 116 || _ent_num == 132) { // Point and Connect Point entities are irrelevant
    } else if (_ent_num == 102) {
      for (int i_curve = 0; i_curve < _parser.read_int(_par[1]); ++i_curve) {
        Read_entity sub_reader(_parser, _parser.read_int(_par[2 + i_curve]), _n_div, _coords);
        Composite_curve sub_curve = sub_reader.read_composite_curve();
        for (auto& c : sub_curve) comp.emplace_back(c.release());
      }
    } else HEXED_THROW("failed to read curve from entity #" + std::to_string(_ent_num));
    return comp;
  }

  // attempts to read any of the entities that derive from `Parametric<2>`
  // Iff `required == true`, throws on failure.
  std::unique_ptr<Parametric<2>> read_surface(bool required = true) const {
    std::unique_ptr<Parametric<2>> ptr;
    merge(ptr, read_plane());
    merge(ptr, read_revolution_surface());
    HEXED_ASSERT(
      !required || ptr,
      "Surface entity #" + std::to_string(_ent_num) + " is not implemented.",
      assert::Not_implemented_error
    );
    return ptr;
  }

  // attempts to read a `Trimmed surface` and returns an empty `optional` on failure
  std::optional<Trimmed_surface> read_trimmed_surface() const {
    if (_ent_num != 144) return {};
    // get surface
    Read_entity surf_reader(_parser, _parser.read_int(_par[1]), _n_div, _coords);
    auto surf = surf_reader.read_surface();
    HEXED_ASSERT(_parser.read_int(_par[2]) == 1, "using outer boundary as boundary curve is not implemented",
                 assert::Not_implemented_error);
    // get trimming curves
    std::vector<Composite_curve> curves;
    for (Int i_curve = 0; i_curve < 1 + _parser.read_int(_par[3]); ++i_curve) {
      Read_entity on_surf(_parser, _parser.read_int(_par[4 + i_curve]), _n_div, _coords);
      HEXED_ASSERT(on_surf._ent_num == 142, "boundary must be a curve on a surface");
      int model_curve = _parser.read_int(on_surf._par[4]);
      if ((_parser.read_int(on_surf._par[5]) == 1 || model_curve == 0) && _parser.read_int(on_surf._par[3]) != 0) {
        HEXED_THROW("parameter-space curves are not implemented", assert::Not_implemented_error);
      } else {
        HEXED_ASSERT(model_curve != 0,
                     "at least one of parameter-space and model-space curve pointers must be defined");
        curves.push_back(Read_entity(_parser, model_curve, _n_div, _coords).read_composite_curve());
      }
    }
    // construct Trimmed_surface
    return {Trimmed_surface(surf.release(), std::move(curves), _n_div)};
  }

  private:
  // assumes that `this` is reading a Transformation Matrix entity and reads it as a coordinate transformation
  Coordinate_change _read_coord() const {
    HEXED_ASSERT(_ent_num == 124, "entity is not a Transformation Matrix");
    Mat<3> translate;
    Mat<3, 3> transform;
    HEXED_ASSERT(_parser.read_int(_dir[6]) == 0, "chaining transformation matrices is not yet implemented");
    auto& par_entry = _parser.entry(Iges_parser::parameter, _parser.read_int(_dir[1]));
    for (int i_dim = 0; i_dim < 3; ++i_dim) {
      for (int j_dim = 0; j_dim < 3; ++j_dim) {
        transform(i_dim, j_dim) = _parser.read_float(par_entry[4*i_dim + j_dim + 1]);
      }
      translate(i_dim) = _unit*_parser.read_float(par_entry[4*i_dim + 3 + 1]);
    }
    return {translate, transform};
  }

  const Iges_parser& _parser;
  const std::vector<std::string>& _dir;
  const std::vector<std::string>& _par;
  Int _ent_num;
  Int _n_div;
  Coordinate_change _coords;
  double _unit;
};

Geom_3d::Geom_3d(std::string file_name, Int n_div) {
  Iges_parser parser(file_name);
  auto dir = parser.section(Iges_parser::directory);
  // read all the trimmed surface entities and ignore everything else
  for (auto& entry : dir) {
    Read_entity read(parser, entry, n_div);
    auto surf = read.read_trimmed_surface();
    if (surf) _surfaces.emplace_back(std::move(*surf));
  }
}

Geom_2d::Geom_2d(std::string file_name, Int n_div) {
  Iges_parser parser(file_name);
  auto dir = parser.section(Iges_parser::directory);
  // read all the curve (not composite curve) entities and ignore everything else
  for (auto& entry : dir) {
    Read_entity read(parser, entry, n_div);
    auto curve = read.read_curve(false);
    if (curve) _curves.emplace_back(curve.release());
  }
}

void Geom_2d::visualize(std::string format, std::string file_name, Int n_div) {
  Int n_nodes = n_div + 1;
  double sz = 1./n_div;
  auto vis = Visualizer::create(format, 3, 1, file_name + "_curves", {}, 0., Visualizer::block);
  for (auto& curve : _curves) {
    Array<double> coords({3, n_nodes});
    for (int i_node = 0; i_node < n_nodes; ++i_node) {
      Mat<3> point = curve->point(Mat<1>{i_node*sz});
      for (int i_dim = 0; i_dim < 3; ++i_dim) coords(i_dim)[i_node] = point(i_dim);
    }
    vis->write_block(coords, Array<double>({0, n_nodes}));
  }
}

double component_bound = default_max_dist/2;

double limit(double dist) {
  return (dist > 0 && dist < default_max_dist) ? dist : default_max_dist;
}

Nearest_point<dyn> recursive_guess_nearest(Mat<> point, double max_distance, double distance_guess,
                                           std::function<Nearest_point<dyn>(Mat<>, double)> bounded_nearest) {
  Nearest_point<dyn> nearest = bounded_nearest(point, distance_guess);
  if ((!nearest.empty() && std::sqrt(nearest.dist_squared()) < distance_guess) || distance_guess >= max_distance) {
    return nearest;
  }
  return recursive_guess_nearest(point, max_distance, distance_guess*2, bounded_nearest);
}

Nearest_point<dyn> guess_nearest(Mat<> point, double max_distance, double distance_guess,
                                 std::function<Nearest_point<dyn>(Mat<> p, double max_dist)> bounded_nearest) {
  // search for the nearest point with `distance_guess` as the maximum distance,
  // and if none are found, recursively double `distance_guess` until `max_distance` is exceeded
  max_distance = limit(max_distance);
  distance_guess = limit(distance_guess);
  for (int i = 0; i < point.size(); ++i) {
    HEXED_ASSERT(-component_bound < point(i) && point(i) < component_bound,
                 format_str(200, "point(%i) == %e is not in bounds", i, point(i)), assert::Numerical_exception);
  }
  return recursive_guess_nearest(point, max_distance, distance_guess, bounded_nearest);
}

Nearest_point<dyn> Geom_2d::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  return guess_nearest(point, max_distance, distance_guess, [this](Mat<> p, double max_dist) {
    Mat<3> p3d = Mat<3>::Zero();
    p3d(Eigen::seqN(0, 2)) = p(Eigen::seqN(0, 2));
    Nearest_point<dyn> nearest(p(Eigen::seqN(0, 2)), max_dist);
    for (auto& curve : _curves) {
      auto param = curve->nearest_params(p3d, [](Mat<1>){return true;}, max_dist);
      if (param.is_feasible) nearest.merge(curve->point(param.params)(Eigen::seqN(0, 2)));
    }
    return nearest;
  });
}

next::Sequence<Mat<3>> Geom_2d::points() {
  // return the endpoints of all curves
  return {
    [this](std::size_t i) {return _curves[i/2]->point(Mat<1>::Constant(i%2));},
    [this]() {return 2*_curves.size();},
  };
}

Nearest_point<dyn> Geom_3d::nearest_point(Mat<> point, double max_distance, double distance_guess) {
  return guess_nearest(point, max_distance, distance_guess, [this](Mat<> p, double max_dist) {
    Nearest_point<dyn> nearest(p, max_dist);
    for (auto& surf : _surfaces) nearest.merge(surf.nearest_point(p, max_dist));
    return nearest;
  });
}

std::vector<double> Geom_3d::intersections(Mat<> start, Mat<> end) {
  Mat<3, 2> endpoints;
  endpoints(all, 0) = start;
  endpoints(all, 1) = end;
  std::vector<double> sects;
  for (auto& surf : _surfaces) {
    std::vector<double> surf_sects = surf.intersections(endpoints);
    sects.insert(sects.end(), surf_sects.begin(), surf_sects.end());
  }
  return sects;
}

next::Sequence<const Tree_curve&> Geom_3d::edges() {
  // concatenate the sequences of bounding curves of all trimmed surfaces
  next::Sequence<const Tree_curve&> e;
  for (auto& surf : _surfaces) e = e + surf.curves();
  return e;
}

void Geom_3d::visualize(std::string format, std::string file_name, Int n_div, bool vis_volume, Mat<3, 2> bounds) {
  Int n_nodes = n_div + 1;
  double sz = 1./n_div;
  {
    auto vis = Visualizer::create(format, 3, 2, file_name + "_surfaces", {"inside"}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      Array<double> discrete({3, n_nodes, n_nodes});
      Array<double> inside({1, n_nodes, n_nodes});
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          Mat<2> params {i*sz, j*sz};
          Mat<3> p = s.surface().point(params);
          for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)(i)[j] = p(i_dim);
          inside(0)(i)[j] = s.is_inside(params);
        }
      }
      vis->write_block(discrete, inside);
    }
  }
  {
    auto vis = Visualizer::create(format, 3, 1, file_name + "_curves", {}, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      for (auto& curve : s.curves()) {
        Array<double> nodes = curve.nodes();
        Array<double> transposed({3, nodes.shape()[0]});
        for (int i = 0; i < nodes.shape()[0]; ++i) {
          for (int i_dim = 0; i_dim < 3; ++i_dim) {
            transposed(i_dim)[i] = nodes(i)[i_dim];
          }
        }
        vis->write_block(transposed, Array<double>({0, n_nodes}));
      }
    }
  }
  if (vis_volume) {
    auto vis = Visualizer::create("default", 3, 3, file_name + "_distance", {"distance"}, 0., Visualizer::block);
    Array<double> coords({3, n_nodes, n_nodes, n_nodes});
    Array<double> dist({1, n_nodes, n_nodes, n_nodes});
    double dist_guess = .1*(bounds(all, 1) - bounds(all, 0)).norm();
    #pragma omp parallel for
    for (int i = 0; i < math::pow(n_nodes, 3); ++i) {
      Mat<3> p;
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        double interp = (i/math::pow(n_nodes, 2 - i_dim)%n_nodes)*sz;
        p(i_dim) = coords(i_dim)[i] = bounds(i_dim, 0) + interp*(bounds(i_dim, 1) - bounds(i_dim, 0));
      }
      Nearest_point np = nearest_point(p, huge, dist_guess);
      dist[i] = (p - np.point()).norm();
    }
    vis->write_block(coords, dist);
  }
}

//! \endcond

}
