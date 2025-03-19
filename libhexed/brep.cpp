#include <optional>
#include <hexed/brep.hpp>
#include <hexed/Iges_parser.hpp>
#include <hexed/Visualizer.hpp>
#include <hexed/utils.hpp>
#include <hexed/constants.hpp>
#include <hexed/Printer.hpp>

namespace hexed::brep {

//! \cond

Coordinate_change::Coordinate_change(Mat<3> translate, Mat<3, 3> transform)
: _translate{translate}, _transform{transform}, _inv{transform.inverse()}
{}

Coordinate_change Coordinate_change::operator()(Coordinate_change that) const {
  return {_translate + _transform*that._translate, _transform*that._transform};
}

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

Mat<2, 1> Line_segment::orig_param_bounds() const {
  HEXED_THROW("not yet implemented for this entity", assert::Not_implemented_error)
  throw;
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

Mat<2, 1> Circular_arc::orig_param_bounds() const {
  HEXED_THROW("not yet implemented for this entity", assert::Not_implemented_error)
  throw;
}

Parametric<2>::Nearest_parameters Plane::nearest_params(Mat<3> p, Constraint is_feasible,
                                                        double max_distance) const {
  Mat<2> params = _vecs.colPivHouseholderQr().solve(p - _origin);
  return {params, is_feasible(params)};
}

std::optional<Mat<2>> Plane::nearest_parameters(Mat<3> p) const {
  Mat<2> params = _vecs.colPivHouseholderQr().solve(p - _origin);
  return {params};
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

Mat<2, 2> Plane::orig_param_bounds() const {
  HEXED_THROW("IGES does not define a parameterization for planes.")
  throw;
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
          double param_diff = segment.radius/std::abs(points_coefs.axial[1]);
          double rsq = points_coefs.radial[0]
                       + (points_coefs.radial[1] + points_coefs.radial[2]*center_param)*center_param;
          double rsq_deriv = points_coefs.radial[1] + 2*points_coefs.radial[2]*center_param;
          double rsq_uncert = (std::abs(rsq_deriv) + points_coefs.radial[2]*param_diff)*param_diff;
          could_intersect = std::abs(radius*radius - rsq) <
                            segment.radius*segment.radius + 2*radius*segment.radius + rsq_uncert;
        } else {
          could_intersect = false;
          std::vector<double> intersection_params;
          for (int layer = 0; layer < 2; ++layer) { // outer, inner
            double target = math::pow(radius - math::sign(layer)*segment.radius, 2);
            double descrim = points_coefs.radial[1]*points_coefs.radial[1]
                             - 4*points_coefs.radial[2]*(points_coefs.radial[0] - target);
            if (descrim > 0) {
              for (int i_sign = 0; i_sign < 2; ++i_sign) {
                intersection_params.insert(intersection_params.begin() + layer + i_sign,
                                           (-points_coefs.radial[1] + math::sign(i_sign)*std::sqrt(descrim))
                                           /(2*points_coefs.radial[2]));
              }
            } else break;
          }
          for (int i_interval = 0; i_interval < int(intersection_params.size()/2); ++i_interval) {
            double bounds [2];
            for (int i = 0; i < 2; ++i) {
              bounds[i] = points_coefs.axial[0] + points_coefs.axial[1]*intersection_params[2*i_interval + i];
            }
            int sign = math::sign(points_coefs.axial[1] > 0);
            if (std::abs(axial - .5*(bounds[0] + bounds[1])) < segment.radius + sign*.5*(bounds[1] - bounds[0])) {
              could_intersect = true;
            }
          }
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

Mat<2, 2> Revolution_surface::orig_param_bounds() const {
  HEXED_THROW("not yet implemented for this entity", assert::Not_implemented_error)
  throw;
}

template class Nurbs<1>;
template class Nurbs<2>;

template <int n_param>
Nurbs<n_param>::Nurbs(std::vector<Array<double>> knots, Array<double> weights, Array<double> control_points,
                      Int n_div, Mat<2, n_param> param_bounds)
: _weights(weights.copy())
, _control_points(control_points.copy())
, _n_div{n_div}
, _orig_bounds{param_bounds}
{
  for (Array<double>& k : knots) _knots.push_back(k.copy());
  HEXED_ASSERT(knots.size() == n_param, "dimensionality mismatch in `knots`")
  HEXED_ASSERT(weights.order() == n_param, "dimensionality mismatch in `weights`")
  HEXED_ASSERT(control_points.order() == n_param + 1, "dimensionality mismatch in `control_points`")
  for (int i_dim = 0; i_dim < n_param; ++i_dim) {
    _n_basis[i_dim] = _weights.shape()[i_dim];
    _degree[i_dim] = _knots[i_dim].size() - _n_basis[i_dim] - 1;
    HEXED_ASSERT(_degree[i_dim] > 0, "degree must be positive")
    HEXED_ASSERT(_control_points.shape()[i_dim] == _n_basis[i_dim],
                 "shape of `weights` and `control_points` don't match")
    _knots[i_dim] -= _orig_bounds(0, i_dim);
    _knots[i_dim] /= _orig_bounds(1, i_dim) - _orig_bounds(0, i_dim);
  }
  HEXED_ASSERT(_control_points.shape()[n_param] == 3, "wrong number of coordinates (should always be 3)")
  double max_sq = 0;
  if constexpr (n_param == 1) {
    #pragma omp parallel for reduction(max:max_sq)
    for (Int i_div = 0; i_div < _n_div; ++i_div) {
      double dist_sq = (point(Mat<1>{double(i_div + 1)/_n_div}) - point(Mat<1>{double(i_div)/_n_div})).squaredNorm()/4;
      max_sq = std::max(max_sq, dist_sq);
    }
  } else {
    #pragma omp parallel for reduction(max:max_sq)
    for (Int i_div = 0; i_div < _n_div; ++i_div) {
      for (Int j_div = 0; j_div < _n_div; ++j_div) {
        Mat<3, 4> verts;
        Mat<3> centroid = Mat<3>::Zero();
        for (int i_vert = 0; i_vert < 2; ++i_vert) {
          for (int j_vert = 0; j_vert < 2; ++j_vert) {
            verts(all, 2*i_vert + j_vert) = point({double(i_div + i_vert)/_n_div, double(j_div + j_vert)/_n_div});
            centroid += verts(all, 2*i_vert + j_vert);
          }
        }
        centroid /= 4;
        for (int i_vert = 0; i_vert < 4; ++i_vert) {
          max_sq = std::max(max_sq, (verts(all, i_vert) - centroid).squaredNorm());
        }
      }
    }
  }
  _max_deriv = std::sqrt(max_sq)*_n_div;
}

template <int n_param>
Mat<3> Nurbs<n_param>::point(Mat<n_param> params) const {
  Int start_knot [n_param];
  std::vector<Array<double>> bases;
  Int n_term = 1;
  for (int i_dim = 0; i_dim < n_param; ++i_dim) {
    bases.push_back(Array<double>({_degree[i_dim] + 1}));
    n_term *= _degree[i_dim] + 1;
    start_knot[i_dim] = _find_knot(i_dim, params[i_dim]);
    Array<double> basis = bases.back()();
    basis = 0;
    basis[0] = 1;
    for (int deg = 1; deg <= _degree[i_dim]; ++deg) {
      for (int j_basis = deg; j_basis >= 0; --j_basis) {
        // note: IGES spec unclear because formulae are for degree k - 1 basis in terms of k - 2
        Array<double> shifted = _knots[i_dim](start_knot[i_dim] - deg, end);
        if (j_basis < deg) {
          basis[j_basis] *= (shifted[j_basis + deg + 1] - params[i_dim])
                            /(shifted[j_basis + deg + 1] - shifted[j_basis + 1]);
        }
        if (j_basis) {
          basis[j_basis] += basis[j_basis - 1]*(params[i_dim] - shifted[j_basis])
                            /(shifted[j_basis + deg] - shifted[j_basis]);
        }
      }
    }
  }
  Mat<3> num = Mat<3>::Zero();
  double denom = 0;
  for (int i_term = 0; i_term < n_term; ++i_term) {
    bool in_bounds = true;
    double basis_fun = 1;
    Int i_cp = 0;
    for (int i_dim = 0; i_dim < n_param; ++i_dim) {
      // this only works for curves and surfaces, no volumes (not that that's a thing anyway)
      static_assert(n_param <= 2);
      int i_basis = (n_param == 2 && !i_dim) ? i_term/(_degree[1] + 1) : i_term%(_degree[i_dim] + 1);
      basis_fun *= bases[i_dim][i_basis];
      Int row = start_knot[i_dim] - _degree[i_dim] + i_basis;
      i_cp += row*_weights.stride(i_dim);
      in_bounds = in_bounds && row < _n_basis[i_dim];
    }
    if (in_bounds) {
      num += _weights[i_cp]*basis_fun*_control_points.reshaped({whatever, 3})(i_cp).vector();
      denom += _weights[i_cp]*basis_fun;
    }
  }
  return num/denom;
}

template <int n_param>
void Nurbs<n_param>::_recursive_nearest(_Nearest_params& nearest, std::array<Int, n_param> start_node, Int size,
                                        Parametric<n_param>::Constraint is_feasible) const {
  constexpr int n_tree = math::pow(2, n_param);
  if (size > 1) {
    double dist [n_tree];
    for (int i_branch = 0; i_branch < n_tree; ++i_branch) {
      Mat<n_param> params;
      for (int i_dim = 0; i_dim < n_param; ++i_dim) {
        Int i_node = start_node[i_dim] + i_branch/math::pow(2, n_param - 1 - i_dim)%2*size/2;
        params(i_dim) = (i_node + size/4.)/_n_div;
      }
      dist[i_branch] = (point(params) - nearest.target).norm();
    }
    for (int i = 0; i < n_tree; ++i) {
      int i_branch = 0;
      for (int j = 1; j < n_tree; ++j) {
        if (dist[j] < dist[i_branch]) i_branch = j;
      }
      std::array<Int, n_param> branch_start = start_node;
      for (int i_dim = 0; i_dim < n_param; ++i_dim) {
        branch_start[i_dim] += i_branch/math::pow(2, n_param - 1 - i_dim)%2*size/2;
      }
      if (math::pow(std::max(0., dist[i_branch] - _max_deriv*size/2/_n_div), 2) < nearest.dist_sq) {
        _recursive_nearest(nearest, branch_start, size/2, is_feasible);
      }
      dist[i_branch] = std::sqrt(huge);
    }
  } else {
    for (int i_vert = 0; i_vert < n_tree; ++i_vert) {
      Mat<n_param> params;
      for (int i_dim = 0; i_dim < n_param; ++i_dim) {
        params(i_dim) = (start_node[i_dim] + i_vert/math::pow(2, n_param - 1 - i_dim)%2)/double(_n_div);
      }
      ++nearest.n_eval;
      if (is_feasible(params)) {
        Mat<3> p = point(params);
        double dist_sq = (p - nearest.target).squaredNorm();
        if (dist_sq < nearest.dist_sq) {
          nearest.params = params;
          nearest.is_feasible = true;
          nearest.dist_sq = dist_sq;
        }
      }
    }
  }
}

template <int n_param>
Parametric<n_param>::Nearest_parameters
Nurbs<n_param>::nearest_params(Mat<3> target, Parametric<n_param>::Constraint is_feasible, double max_distance) const {
  _Nearest_params nearest {target, Mat<n_param>::Zero(), false, max_distance*max_distance, 0};
  std::array<Int, n_param> start_node;
  start_node.fill(0);
  _recursive_nearest(nearest, start_node, _n_div, is_feasible);
  return {nearest.params, nearest.is_feasible};
}

template <int n_param>
Int Nurbs<n_param>::_find_knot(int i_dim, double param) const {
  // O(log n) binary search
  Int low = _degree[i_dim];
  Int high = _knots[i_dim].size() - 1 - _degree[i_dim];
  while (high - low > 1) {
    Int diff = (high - low)/2;
    if (_knots[i_dim][low + diff] < param) {
      low += diff;
    } else if (_knots[i_dim][high - diff] > param) {
      high -= diff;
    // if `param` is exactly equal to a knot, the following cases will be selected
    } else if (_knots[i_dim][low + diff] < .5*(_knots[i_dim][high] + _knots[i_dim][low])) {
      low += diff;
    } else {
      high -= diff;
    }
  }
  return low;
}

template <int n_param>
std::vector<typename Parametric<n_param>::Intersection_parameters>
Nurbs<n_param>::intersection_params(Mat<3, 2> points) const {
  return {};
}

template <int n_param>
Mat<2, n_param> Nurbs<n_param>::orig_param_bounds() const {
  return _orig_bounds;
}

Trimmed_surface::Trimmed_surface(Parametric<2>* surface, std::vector<Composite_curve>&& curves,
                                 std::vector<bool> is_model_space, Int n_div)
: _n_div{n_div}
, _sz{1./_n_div}
, _surf{surface}
, _levels{math::log(2, _n_div)}
, _excession({_levels + 1, 2})
, _excession_epsilon({2})
{
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      Array<double> boundary({_n_div + 1, 3});
      #pragma omp parallel for
      for (Int i_node = 0; i_node < _n_div + 1; ++i_node) {
        Mat<2> params;
        params(i_dim) = sign;
        params(!i_dim) = i_node*_sz;
        boundary(i_node).vector() = _surf->point(params);
      }
      _extremal_boundaries.emplace_back(boundary.copy());
    }
  }
  Mat<3> ref_point = _surf->point(Mat<2>{.5, .5});
  double ref_dist = 0;
  #pragma omp parallel for reduction(max:ref_dist)
  for (Int i_node = 0; i_node < _n_div + 1; ++i_node) {
    for (Int j_node = 0; j_node < _n_div + 1; ++j_node) {
      ref_dist = std::max(ref_dist, (_surf->point(Mat<2>{i_node*_sz, j_node*_sz}) - ref_point).norm());
    }
  }
  _excession_epsilon = Array<double>::make(1e-3*ref_dist*_sz, 1e-3*_n_div);
  _excession = 0; // this will set `_excession(_levels)`, which is not set in the following loop
  for (int level = 0; level < _levels; ++level) {
    Int n_outer = math::pow(2, level);
    Int n_inner = _n_div/n_outer;
    for (Int i_outer = 0; i_outer < n_outer; ++i_outer) {
      for (Int j_outer = 0; j_outer < n_outer; ++j_outer) {
        Array<double> dist_nrml({4, 2, 3});
        Array<double> average({2, 3});
        average = 0;
        for (int i_vert = 0; i_vert < 4; ++i_vert) {
          Mat<2> params {(i_outer + i_vert/2*n_inner)*_sz, (j_outer + i_vert%2*n_inner)*_sz};
          dist_nrml(i_vert)(0).vector() = _surf->point(params);
          dist_nrml(i_vert)(1).vector() = normal(params);
          average += dist_nrml(i_vert)/4.;
        }
        double approx_radius [2] {};
        for (int i_vert = 0; i_vert < 4; ++i_vert) {
          for (int i = 0; i < 2; ++i) {
            approx_radius[i] = std::max(approx_radius[i], (dist_nrml(i_vert)(i) - average(i)).vector().norm());
          }
        }
        for (Int i_inner = 0; i_inner < n_inner + 1; ++i_inner) {
          for (Int j_inner = 0; j_inner < n_inner + 1; ++j_inner) {
            Mat<2> params {(i_outer + i_inner)*_sz, (j_outer + j_inner)*_sz};
            Array<double> dn({2, 3});
            dn(0).vector() = _surf->point(params);
            dn(1).vector() = normal(params);
            for (int i = 0; i < 2; ++i) {
              double r = (dn(i) - average(i)).vector().norm();
              _excession(level)[i] = std::max(_excession(level)[i],
                                              (r - approx_radius[i])/(approx_radius[i] + _excession_epsilon[i]));
            }
          }
        }
      }
    }
  }
  // discretize curves into polygonal segments in parameter space
  std::vector<std::vector<std::vector<Mat<2>>>> discrete_curves;
  for (int i_composite = 0; i_composite < int(curves.size()); ++i_composite) {
    auto& composite = curves[i_composite];
    discrete_curves.emplace_back();
    auto& disc_curve = discrete_curves.back();
    for (auto& curve : composite) {
      disc_curve.emplace_back();
      auto& param_nodes = disc_curve.back();
      Array<double> phys_nodes({_n_div + 1, 3});
      double mean_squared_dist = 0;
      Mat<3> start;
      if (is_model_space[i_composite]) {
        start = curve->point(Mat<1>{0.});
        for (Int i_node = 0; i_node < _n_div + 1; ++i_node) {
          Mat<3> pt = curve->point(Mat<1>{i_node*_sz});
          mean_squared_dist += (pt - start).squaredNorm();
          phys_nodes(i_node).vector() = pt;
          Mat<2> params = _nearest_params(pt);
          param_nodes.push_back(params);
        }
      } else {
        Mat<2, 2> op = _surf->orig_param_bounds();
        auto get_params = [&curve, op](double t) {
          Mat<2> p = curve->point(Mat<1>{t})(Eigen::seqN(0, 2));
          for (int i_dim = 0; i_dim < 2; ++i_dim) {
            p(i_dim) = (p(i_dim) - op(0, i_dim))/(op(1, i_dim) - op(0, i_dim));
          }
          return p;
        };
        start = _surf->point(get_params(0.));
        for (Int i_node = 0; i_node < _n_div + 1; ++i_node) {
          Mat<2> params = get_params(i_node*_sz);
          Mat<3> pt = _surf->point(params);
          mean_squared_dist += (pt - start).squaredNorm();
          phys_nodes(i_node).vector() = pt;
          param_nodes.push_back(params);
        }
      }
      mean_squared_dist /= n_div + 1;
      if ((curve->point(Mat<1>{1.}) - start).squaredNorm() < .1*mean_squared_dist) {
        _curves.emplace_back(phys_nodes(0, n_div/2 + 1).copy(), 4);
        _curves.emplace_back(phys_nodes(n_div/2, n_div + 1).copy(), 4);
      } else {
        _curves.emplace_back(phys_nodes.copy(), 4);
      }
    }
    int reversal = 0;
    double continuity_error = huge;
    int n_curves = disc_curve.size();
    for (int test_reversal = 0; test_reversal < math::pow<int>(2, n_curves); ++test_reversal) {
      Array<double> endpoints({n_curves, 2, 3});
      for (int i_curve = 0; i_curve < n_curves; ++i_curve) {
        bool reverse = test_reversal%math::pow(2, i_curve + 1)/math::pow(2, i_curve);
        for (int i_end = 0; i_end < 2; ++i_end) {
          Int i_point = (i_end + reverse)%2*(disc_curve[i_curve].size() - 1);
          endpoints(i_curve)(i_end).vector() = _surf->point(disc_curve[i_curve][i_point]);
        }
      }
      double err = 0;
      for (int i_curve = 0; i_curve < n_curves; ++i_curve) {
        err += (endpoints((i_curve + 1)%n_curves)(0) - endpoints(i_curve)(1)).vector().norm();
      }
      if (err < continuity_error) {
        continuity_error = err;
        reversal = test_reversal;
      }
    }
    for (int i_curve = 0; i_curve < n_curves; ++i_curve) {
      if (reversal%math::pow(2, i_curve + 1)/math::pow(2, i_curve)) {
        std::reverse(disc_curve[i_curve].begin(), disc_curve[i_curve].end());
      }
    }
  }
  // initialize parameter-space curves with discretiation
  _initialize(discrete_curves);
}

void Trimmed_surface::_initialize(std::vector<std::vector<std::vector<Mat<2>>>>& curves) {
  // compute bounds of discrete nodes in parameter space
  Mat<2, 2> bounds;
  bounds << huge, -huge, huge, -huge;
  bool set = false;
  for (auto& loop : curves) {
    for (auto& curve : loop) {
      set = set || !curve.empty();
      for (Mat<2> node : curve) {
        bounds(all, 0) = bounds(all, 0).cwiseMin(node);
        bounds(all, 1) = bounds(all, 1).cwiseMax(node);
      }
    }
  }
  if (set) bounds(all, 1) = bounds(all, 1).cwiseMax(bounds(all, 0) + Mat<2>{_sz, _sz});
  else bounds << 0, 1, 0, 1;
  // reparameterize surface to contain bounds
  bounds = _surf->reparameterize(bounds);
  for (int i_direction = 0; i_direction < 3; ++i_direction) _param_segments[i_direction].resize(_n_div);
  for (auto& loop : curves) {
    std::vector<Mat<2>> all_nodes;
    std::vector<std::array<Int, 2>> curve_endpoints;
    Int last_endpoint = 0;
    for (auto& nodes : loop) {
      curve_endpoints.push_back({last_endpoint, last_endpoint += nodes.size()});
      all_nodes.insert(all_nodes.end(), nodes.begin(), nodes.end());
    }
    Int n_nodes = all_nodes.size();
    // apply reparameterization to nodes
    for (Int i_node = 0; i_node < n_nodes; ++i_node) {
      all_nodes[i_node] = (all_nodes[i_node] - bounds(all, 0)).cwiseQuotient(bounds(all, 1) - bounds(all, 0));
    }
    // correct periodic seam errors
    for (Int i_node = n_nodes, changed = false; (i_node < 2*n_nodes) || changed; ++i_node, changed = false) {
      for (int i_dim = 0; i_dim < 2; ++i_dim) {
        for (int sign : {-1, 1}) {
          if (std::abs(all_nodes[i_node%n_nodes](i_dim) + sign - all_nodes[(i_node - 1)%n_nodes](i_dim)) <
              std::abs(all_nodes[i_node%n_nodes](i_dim)        - all_nodes[(i_node - 1)%n_nodes](i_dim))) {
            all_nodes[i_node%n_nodes](i_dim) += sign;
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
      for (Mat<2> node : all_nodes) {
        n_less += node(i_dim) < -_sz;
        n_greater += node(i_dim) > 1 + _sz;
      }
      int sign = (n_less > n_nodes/2) - (n_greater > n_nodes/2);
      for (Mat<2>& node : all_nodes) node(i_dim) += sign;
    }
    // check if there's any two curves that are exactly the same, but reversed,
    // and also on one of the parameter extremes,
    // then that's a clue that they should actually be on opposite sides but aren't
    // because there is both a periodic seam and a pole.
    int n_curve = loop.size();
    for (int i_curve = 0; i_curve < n_curve; ++i_curve) {
      Int size0 = curve_endpoints[i_curve][1] - curve_endpoints[i_curve][0];
      for (int j_curve = 0; j_curve < i_curve; ++j_curve) {
        Int size1 = curve_endpoints[j_curve][1] - curve_endpoints[j_curve][0];
        if (size1 == size0) {
          for (int i_dim = 0; i_dim < 2; ++i_dim) {
            Int n_duplicate [2] {};
            for (Int i_node = 0; i_node < size0; ++i_node) {
              for (int side : {0, 1}) {
                Int node0 = curve_endpoints[i_curve][0] + i_node;
                Int node1 = curve_endpoints[j_curve][0] + size0 - 1 - i_node;
                n_duplicate[side] += std::abs(all_nodes[node0](i_dim) - side) < _sz &&
                                     std::abs(all_nodes[node1](i_dim) - side) < _sz &&
                                     std::abs(all_nodes[node0](!i_dim) - all_nodes[node1](!i_dim)) < _sz;
              }
            }
            for (int side : {0, 1}) {
              // at least 2 of the nodes might not match up because of singularities,
              // so let's say 4 can differ to be safe
              if (n_duplicate[side] > size0 - 4) {
                for (Int i_node = curve_endpoints[i_curve][0]; i_node < curve_endpoints[i_curve][1]; ++i_node) {
                  all_nodes[i_node](i_dim) -= math::sign(side);
                }
              }
            }
          }
        }
      }
    }
    if (!all_nodes.empty()) {
      all_nodes.insert(all_nodes.end(), all_nodes.front());
      ++n_nodes;
    }
    // compute parametric segments
    for (int i_direction = 0; i_direction < 3; ++i_direction) {
      Mat<2, 2> transform = _transform_mat(i_direction);
      std::vector<Int> abscissa;
      std::vector<double> ordinate;
      Mat<2> prev_params {-1., 0.};
      Int prev_absc = -1;
      if (n_nodes) for (Int i_node = 0; i_node <= n_nodes; ++i_node) {
        Mat<2> params;
        if (i_node == n_nodes && !abscissa.empty()) params << abscissa.front(), ordinate.front();
        else {
          params = transform*all_nodes[i_node];
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
        HEXED_ASSERT(abscissa.front() == abscissa.back(),
                     format_str(200, "parametric representation is not closed (%li != %li)",
                                abscissa.front(), abscissa.back()))
        for (Int i_segment = 0; i_segment < Int(abscissa.size()) - 1; ++i_segment) {
          HEXED_ASSERT(std::abs(abscissa[i_segment] - abscissa[i_segment + 1]) == 1, "invalid step size");
          bool reverse = abscissa[i_segment] > abscissa[i_segment + 1];
          HEXED_ASSERT(0 <= abscissa[i_segment + reverse] && abscissa[i_segment + reverse] < _n_div,
                       format_str(200, "segment index %i out of bounds", abscissa[i_segment + reverse]));
          _param_segments[i_direction][abscissa[i_segment + reverse]].push_back({ordinate[i_segment +  reverse],
                                                                                 ordinate[i_segment + !reverse]});
        }
      }
    }
  }
}

Mat<2, 2> Trimmed_surface::_transform_mat(int i_direction) const {
  Mat<2> dir;
  switch (i_direction) {
    case 0: dir << 1., 0.; break;
    case 1: dir << std::sqrt(3.)/2, -.5; break;
    case 2: dir << .5, -std::sqrt(3.)/2; break;
    default: HEXED_THROW("`i_direction` must be in [0, 3)")
  }
  dir /= dir(0) - dir(1);
  Mat<2, 2> trans;
  trans <<
    dir(0), -dir(1),
    dir(1),  dir(0);
  return trans;
}

bool Trimmed_surface::is_inside(Mat<2> params) const {
  // fudge it a little so the outer boundary is always inside
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    double tol = 1e-3/_n_div;
    if (params(i_dim) < -tol || params(i_dim) > 1 + tol) return false;
    params(i_dim) = std::max(tol, std::min(1 - tol, params(i_dim)));
  }
  Int n_intersections = 0;
  double diff = huge;
  for (int i_direction = 0; i_direction < 3; ++i_direction) {
    Mat<2> p = _transform_mat(i_direction)*params;
    // count the number of segmements intersected by a ray in the positive `p(1)` direction
    Int i_seg = floor(p(0)*_n_div);
    double max_diff = 0;
    for (Int j_seg : {i_seg - 1, i_seg, i_seg + 1}) if (0 <= j_seg && j_seg < _n_div) {
      for (Mat<2> seg : _param_segments[i_direction][j_seg]) {
        if (std::min(seg(0), seg(1)) <= p(1) && p(1) <= std::max(seg(0), seg(1))) {
          max_diff = std::max(max_diff, std::abs(seg(0) - seg(1)));
        }
      }
    }
    if (max_diff < diff) {
      diff = max_diff;
      n_intersections = 0;
      for (Mat<2> seg : _param_segments[i_direction][i_seg]) {
        n_intersections += p(1) < seg(0) + (p(0)*_n_div - i_seg)*(seg(1) - seg(0));
      }
    }
  }
  // the point is inside iff the number of intesections is odd
  return n_intersections%2;
}

next::Sequence<const Tree_curve&> Trimmed_surface::curves() const {
  return next::Sequence<const Tree_curve&>::vector_view(_curves);
}

void Trimmed_surface::_recursive_nearest(Nearest_point<3>& nearest, Mat<2>& best_params, Int i_start, Int j_start,
                                         Int level, bool check_inside) const {
  Int n_panel = _n_div/math::pow(2, level);
  Mat<3> point = nearest.reference();
  Array<double> dist_nrml({2, 2, 2, 3});
  Array<double> average({2, 3});
  average = 0;
  for (int i_vertex = 0; i_vertex < 2; ++i_vertex) {
    for (int j_vertex = 0; j_vertex < 2; ++j_vertex) {
      Mat<2> params {(i_start + i_vertex*n_panel)*_sz, (j_start + j_vertex*n_panel)*_sz};
      Mat<3> n = normal(params);
      Mat<3> d = _surf->point(params) - point;
      dist_nrml(i_vertex)(j_vertex)(0).vector() = d;
      dist_nrml(i_vertex)(j_vertex)(1).vector() = n;
      average += dist_nrml(i_vertex)(j_vertex);
    }
  }
  average /= 4;
  double radii [2] {};
  for (int i_vertex = 0; i_vertex < 2; ++i_vertex) {
    for (int j_vertex = 0; j_vertex < 2; ++j_vertex) {
      for (int i = 0; i < 2; ++i) {
        radii[i] = std::max(radii[i], (dist_nrml(i_vertex)(j_vertex)(i) - average(i)).vector().norm());
      }
    }
  }
  for (int i = 0; i < 2; ++i) radii[i] += _excession(level)[i]*(radii[i] + _excession_epsilon[i]);
  double norm = average(1).vector().norm() + 1e-12;
  Mat<3> unit_avg = average(1).vector()/norm;
  bool compute;
  if (average(0).vector().norm() > radii[0] + std::sqrt(nearest.dist_squared())) {
    compute = false;
  } else if (radii[1] > norm - 1e-3) {
    compute = true;
  } else {
    double scale = average(0).vector().dot(unit_avg)/norm;
    average(1) *= scale;
    double sin = radii[1]/norm;
    double cos = std::sqrt(1 - sin*sin);
    average(1).vector() += (average(1) - average(0)).vector().norm()*sin/cos*math::sign(scale > 0)*unit_avg;
    double r = radii[0] + (radii[1] + 1e-3)*average(1).vector().norm()/norm;
    compute = (average(0) - average(1)).vector().norm() < r;
  }
  if (compute) {
    if (n_panel == 1) {
      Mat<2> params;
      for (int i_triangle = 0; i_triangle < 2; ++i_triangle) {
        Mat<3, 2> lhs;
        lhs(all, 0) = (dist_nrml(!i_triangle)(i_triangle)(0) - dist_nrml(i_triangle)(i_triangle)(0)).vector();
        lhs(all, 1) = (dist_nrml(i_triangle)(!i_triangle)(0) - dist_nrml(i_triangle)(i_triangle)(0)).vector();
        Mat<2> soln = lhs.householderQr().solve(-dist_nrml(i_triangle)(i_triangle)(0).vector());
        if (!(soln(0) >= 0 && soln(1) >= 0 && soln(0) + soln(1) <= 1)) {
          double dist_sq = huge;
          for (int i_edge = 0; i_edge < 3; ++i_edge) {
            Mat<3> start = dist_nrml(i_triangle != (i_edge == 1))(i_triangle != (i_edge == 2))(0).vector();
            int i_end = (i_edge + 1)%3;
            Mat<3> diff  = dist_nrml(i_triangle != (i_end == 1))(i_triangle != (i_end == 2))(0).vector() - start;
            double interp = std::max(0., std::min(1., -start.dot(diff)/diff.squaredNorm()));
            double d = (start + interp*diff).squaredNorm();
            if (d < dist_sq) {
              dist_sq = d;
              switch (i_edge) {
                case 0: soln(0) = interp; soln(1) = 0; break;
                case 1: soln(0) = 1 - interp; soln(1) = interp; break;
                case 2: soln(0) = 0; soln(1) = 1 - interp; break;
              }
            }
          }
        }
        params << (i_start + i_triangle)*_sz, (j_start + i_triangle)*_sz;
        params -= math::sign(i_triangle)*soln*_sz;
        if (!check_inside || is_inside(params)) {
          if (nearest.merge(point + dist_nrml(i_triangle)(i_triangle)(0).vector() + lhs*soln)) best_params = params;
        }
      }
    } else {
      for (int i_subsect = 0; i_subsect < 2; ++i_subsect) {
        for (int j_subsect = 0; j_subsect < 2; ++j_subsect) {
          _recursive_nearest(nearest, best_params, i_start + i_subsect*n_panel/2, j_start + j_subsect*n_panel/2,
                             level + 1, check_inside);
        }
      }
    }
  }
}

Nearest_point<3> Trimmed_surface::nearest_point(Mat<3> point, double max_dist) const {
  Nearest_point<3> nearest(point, max_dist);
  Mat<2> best_params; // unused
  // first check all local nearest points in the interior of the surface
  _recursive_nearest(nearest, best_params, 0, 0, 0, true);
  // then check the nearest point on all the boundary curves
  for (auto& curve : _curves) {
    auto index = curve.nearest_point(point, 1.01*std::sqrt(nearest.dist_squared()));
    if (index.index > -1) nearest.merge(curve.nodes()(index.index).vector());
  }
  return nearest;
}

Mat<2> Trimmed_surface::_nearest_params(Mat<3> point) const {
  Nearest_point<3> nearest(point, default_max_dist);
  Mat<2> best_params = Mat<2>::Zero();
  auto par = _surf->nearest_parameters(point);
  if (par.has_value()) {
    best_params = par.value();
    nearest.merge(_surf->point(par.value()));
  } else {
    _recursive_nearest(nearest, best_params, 0, 0, 0, false);
  }
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign = 0; sign < 2; ++sign) {
      auto& curve = _extremal_boundaries[2*i_dim + sign];
      auto index = curve.nearest_point(point, 1.01*std::sqrt(nearest.dist_squared()));
      if (index.index > -1) {
        if (nearest.merge(curve.nodes()(index.index).vector())) {
          best_params(i_dim) = sign;
          best_params(!i_dim) = index.interp_index*_sz;
        }
      }
    }
  }
  return best_params;
}

std::vector<double> Trimmed_surface::intersections(Mat<3, 2> endpoints) const {
  std::vector<double> sects;
  auto sect_params = _surf->intersection_params(endpoints);
  for (auto params : sect_params) {
    if (is_inside(params.params)) sects.push_back(params.interp_coef);
  }
  return sects;
}

Mat<3> Trimmed_surface::normal(Mat<2> params) const {
  Mat<3, 2> diffs = Mat<3, 2>::Zero();
  for (int i_dim = 0; i_dim < 2; ++i_dim) {
    for (int sign : {-1, 1}) {
      Mat<2> p = params;
      p(i_dim) += sign*.5*_sz;
      Mat<3> diff_point;
      if (p(i_dim) < 0 || p(i_dim) > 1) {
        Mat<3> point0 = _surf->point(params);
        Mat<3> point1 = _surf->point(1.5*params - .5*p);
        Mat<3> point2 = _surf->point(2*params - p);
        diff_point = 6*point0 - 8*point1 + 3*point2;
      } else {
        diff_point = _surf->point(p);
      }
      diffs(all, i_dim) += sign*diff_point;
    }
  }
  return diffs(all, 0).cross(diffs(all, 1)).normalized();
}

Mat<3> Trimmed_surface::point(Mat<2> params) const {
  return _surf->point(params);
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

  Read_entity make_reader(Int line, Coordinate_change change_to = {}) const {
    return {_parser, line, _n_div, _coords(change_to)};
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

  std::unique_ptr<Line_segment> read_line_segment() const {
    if (_ent_num != 110) return {};
    Mat<3, 2> endpts;
    for (int i = 0; i < 6; ++i) endpts(i) = _unit*_parser.read_float(_par[1 + i]);
    for (int i = 0; i < 2; ++i) endpts(all, i) = _coords.to_model(endpts(all, i));
    return std::make_unique<Line_segment>(endpts);
  }

  std::unique_ptr<Parametric<1>> read_circular_arc() const {
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
    auto axis = make_reader(_parser.read_int(_par[1])).read_line_segment();
    auto generatrix = make_reader(_parser.read_int(_par[2])).read_curve();
    return std::make_unique<Revolution_surface>(generatrix.release(), *axis, _n_div,
                                                _parser.read_float(_par[3]), _parser.read_float(_par[4]));
  }

  std::unique_ptr<Nurbs<1>> read_nurbs_curve() const {
    if (_ent_num != 126) return {};
    Int n_basis = _parser.read_int(_par[1]) + 1;
    Int degree = _parser.read_int(_par[2]);
    Int i = 7;
    std::vector<Array<double>> knots;
    knots.emplace_back(std::vector<Int>{1 + n_basis + degree});
    for (int i_knot = 0; i_knot < knots[0].size(); ++i_knot) {
      knots[0][i_knot] = _parser.read_float(_par[i++]);
    }
    Array<double> weights({n_basis});
    for (int i_weight = 0; i_weight < n_basis; ++i_weight) {
      weights[i_weight] = _parser.read_float(_par[i++]); // note transposed
    }
    Array<double> control_points({n_basis, 3});
    for (int i_point = 0; i_point < n_basis; ++i_point) {
      for (int i_dim = 0; i_dim < 3; ++i_dim) {
        control_points(i_point)[i_dim] = _unit*_parser.read_float(_par[i++]);
      }
      control_points(i_point).vector() = _coords.to_model(control_points(i_point).vector());
    }
    Mat<2, 1> param_bounds;
    param_bounds << _parser.read_float(_par[i]), _parser.read_float(_par[i + 1]);
    return std::make_unique<Nurbs<1>>(std::move(knots), weights(), control_points(), _n_div, param_bounds);
  }

  std::unique_ptr<Nurbs<2>> read_nurbs_surface() const {
    if (_ent_num != 128) return {};
    std::vector<Int> n_basis {_parser.read_int(_par[1]) + 1, _parser.read_int(_par[2]) + 1};
    std::vector<Int> degree  {_parser.read_int(_par[3]), _parser.read_int(_par[4])};
    Int i = 10;
    std::vector<Array<double>> knots;
    for (int i_dim = 0; i_dim < 2; ++i_dim) {
      knots.emplace_back(std::vector<Int>{1 + n_basis[i_dim] + degree[i_dim]});
      for (int i_knot = 0; i_knot < knots[i_dim].size(); ++i_knot) {
        knots[i_dim][i_knot] = _parser.read_float(_par[i++]);
      }
    }
    Array<double> weights(n_basis);
    for (int i_weight = 0; i_weight < n_basis[1]; ++i_weight) {
      for (int j_weight = 0; j_weight < n_basis[0]; ++j_weight) {
        weights(j_weight)[i_weight] = _parser.read_float(_par[i++]); // note transposed
      }
    }
    Array<double> control_points({n_basis[0], n_basis[1], 3});
    for (int i_point = 0; i_point < n_basis[1]; ++i_point) {
      for (int j_point = 0; j_point < n_basis[0]; ++j_point) {
        for (int i_dim = 0; i_dim < 3; ++i_dim) {
          control_points(j_point)(i_point)[i_dim] = _unit*_parser.read_float(_par[i++]);
        }
        control_points(j_point)(i_point).vector() = _coords.to_model(control_points(j_point)(i_point).vector());
      }
    }
    Mat<2, 2> param_bounds;
    param_bounds <<
      _parser.read_float(_par[i + 0]), _parser.read_float(_par[i + 2]),
      _parser.read_float(_par[i + 1]), _parser.read_float(_par[i + 3]);
    return std::make_unique<Nurbs<2>>(std::move(knots), weights(), control_points(), _n_div, param_bounds);
  }

  // Attempts to read any of the entities that derive from `Parametric<1>`.
  // Iff `required == true`, throws on failure.
  std::unique_ptr<Parametric<1>> read_curve(bool required = true) const {
    std::unique_ptr<Parametric<1>> ptr;
    merge(ptr, read_line_segment());
    merge(ptr, read_circular_arc());
    merge(ptr, read_nurbs_curve());
    HEXED_ASSERT(
      !required || ptr,
      "Curve entity #" + std::to_string(_ent_num) + " is not implemented.",
      assert::Not_implemented_error
    );
    return ptr;
  }

  // attempts to read a composite curve (collection of curves that share endpoints) and throws on failure
  Composite_curve read_composite_curve() const {
    auto ptr = read_curve(false);
    Composite_curve comp;
    if (ptr) {
      comp.emplace_back(ptr.release());
    } else if (_ent_num == 116 || _ent_num == 132) { // Point and Connect Point entities are irrelevant
    } else if (_ent_num == 102) {
      for (int i_curve = 0; i_curve < _parser.read_int(_par[1]); ++i_curve) {
        Composite_curve sub_curve = make_reader(_parser.read_int(_par[2 + i_curve])).read_composite_curve();
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
    merge(ptr, read_nurbs_surface());
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
    auto surf = make_reader(_parser.read_int(_par[1])).read_surface();
    HEXED_ASSERT(_parser.read_int(_par[2]) == 1, "using outer boundary as boundary curve is not implemented",
                 assert::Not_implemented_error);
    // get trimming curves
    std::vector<Composite_curve> curves;
    std::vector<bool> model_space;
    for (Int i_curve = 0; i_curve < 1 + _parser.read_int(_par[3]); ++i_curve) {
      Read_entity on_surf = make_reader(_parser.read_int(_par[4 + i_curve]));
      HEXED_ASSERT(on_surf._ent_num == 142, "boundary must be a curve on a surface");
      int model_curve = _parser.read_int(on_surf._par[4]);
      int param_curve = _parser.read_int(on_surf._par[3]);
      if ((_parser.read_int(on_surf._par[5]) == 1 || model_curve == 0) && param_curve != 0) {
        model_space.push_back(false);
        Coordinate_change unscale(Mat<3>::Zero(), Mat<3, 3>::Identity()/_unit);
        curves.push_back(make_reader(param_curve, unscale).read_composite_curve());
      } else {
        HEXED_ASSERT(model_curve != 0,
                     "at least one of parameter-space and model-space curve pointers must be defined");
        model_space.push_back(true);
        curves.push_back(make_reader(model_curve).read_composite_curve());
      }
    }
    // construct Trimmed_surface
    return {Trimmed_surface(surf.release(), std::move(curves), std::move(model_space), _n_div)};
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

next::Sequence<const Trimmed_surface&> Geom_3d::surfaces() {
  return next::Sequence<const Trimmed_surface&>::vector_view(_surfaces);
}

void Geom_3d::visualize(std::string format, std::string file_name, Int n_div, bool vis_volume, Mat<3, 2> bounds) {
  Int n_nodes = n_div + 1;
  double sz = 1./n_div;
  {
    std::vector<std::string> var_names {"normal0", "normal1", "normal2", "inside"};
    auto vis = Visualizer::create(format, 3, 2, file_name + "_surfaces", var_names, 0., Visualizer::block);
    for (auto& s : _surfaces) {
      Array<double> discrete({3, n_nodes, n_nodes});
      Array<double> data({4, n_nodes, n_nodes});
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          Mat<2> params {i*sz, j*sz};
          Mat<3> p = s.surface().point(params);
          for (int i_dim = 0; i_dim < 3; ++i_dim) discrete(i_dim)(i)[j] = p(i_dim);
          Mat<3> n = s.normal(params);
          for (int i_dim = 0; i_dim < 3; ++i_dim) data(i_dim)(i)[j] = n(i_dim);
          data(3)(i)[j] = s.is_inside(params);
        }
      }
      vis->write_block(discrete, data);
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
