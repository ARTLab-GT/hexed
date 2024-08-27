#ifndef HEXED_BREP_HPP_
#define HEXED_BREP_HPP_

#include <memory>
#include <functional>
#include "math.hpp"
#include "constants.hpp"
#include "Tree_curve.hpp"
#include "Surface_geom.hpp"

namespace hexed::brep {

constexpr double default_max_dist = std::sqrt(huge);

template <int n_param>
class Parametric {
  public:
  struct Nearest_parameters {
    Mat<n_param> params;
    bool is_feasible;
  };
  typedef std::function<bool(Mat<n_param>)> Constraint;
  virtual Mat<3> point(Mat<n_param> params) const = 0;
  virtual Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                            double max_distance) const = 0;
  Mat<3> nearest_point(Mat<3> p) const {
    return point(nearest_params(p, [](Mat<n_param>){return true;}, default_max_dist).params);
  }
  virtual void reparameterize(Mat<n_param, 2> bounds) {}
};

class Line_segment : public Parametric<1> {
  public:
  inline Line_segment(Mat<3, 2> endpoints) : _endpoints{endpoints} {}
  inline Mat<3> point(Mat<1> params) const override {return _endpoints*Mat<2>{1. - params(0), params(0)};}
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                    double max_distance) const override;
  private:
  Mat<3, 2> _endpoints;
};

class Circular_arc : public Parametric<1> {
  public:
  Circular_arc(Mat<3> center, double radius, double start_angle, double end_angle);
  Mat<3> point(Mat<1> params) const override;
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                    double max_distance) const override;
  private:
  Mat<3> _center;
  double _radius;
  double _start_angle;
  double _end_angle;
};

class Plane : public Parametric<2> {
  public:
  inline Plane(Mat<3> origin, Mat<3, 2> coord_vectors) : _origin{origin}, _vecs{coord_vectors} {}
  void reparameterize(Mat<2, 2> bounds) override {
    _origin = point(bounds(all, 0));
    _vecs = _vecs*(bounds(all, 1) - bounds(all, 0)).asDiagonal();
  }
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                    double max_distance) const override;
  inline Mat<3> point(Mat<2> params) const override {return _origin + _vecs*params;}
  private:
  Mat<3> _origin;
  Mat<3, 2> _vecs;
};

#if 0
struct Trans_mat {
  Mat<3, 3> transform = Mat<3, 3>::Identity();
  Mat<3> translate = Mat<3>::Zero();
};

template <int n_param>
class Entity {
  public:
  struct Nearest_params {
    Mat<n_param> params;
    bool is_feasible;
  };
  typedef std::function<bool(Mat<n_param>)> Constraint;
  virtual ~Entity() = default;
  Nearest_params nearest_params(Mat<3> p, Constraint is_feasible, double max_distance = default_max_dist) const {
    return temp_nearest_params(inv_convert(p), is_feasible, max_distance);
  }
  virtual Mat<3> nearest_point(Mat<3> p, double max_distance = default_max_dist) const {
    Nearest_params params = nearest_params(p, [](Mat<n_param>){return true;}, max_distance);
    HEXED_ASSERT(params.is_feasible, "no feasible nearest point found");
    return point(params.params);
  }
  Mat<3> point(Mat<n_param> params) const {
    return _convert(temp_point(params));
  }
  virtual void reparameterize(Mat<n_param, 2> bounds) {}
  Int n_div = math::pow(2, 13);
  double scale = 1.;
  constexpr static double default_max_dist = std::sqrt(huge);
  Trans_mat trans_mat;
  Mat<3> _convert(Mat<3> p) const {
    return scale*(trans_mat.translate + trans_mat.transform*p);
  }
  Mat<3> inv_convert(Mat<3> p) const {
    return trans_mat.transform.colPivHouseholderQr().solve(p/scale - trans_mat.translate);
  }
  virtual Nearest_params temp_nearest_params(Mat<3>, Constraint is_feasible, double max_distance) const {
    return {Mat<n_param>::Zero(), false};
  };
  virtual Mat<3> temp_point(Mat<n_param>) const = 0;
  private:
};

class Revolution_surface : public Entity<2> {
  public:
  Revolution_surface(Entity<1>*, Line_segment, double start_angle = 0, double end_angle = 2*constants::pi);
  Mat<3> rotate(Mat<3>, double angle) const;
  std::unique_ptr<Entity<1>> generatrix;
  Line_segment axis;
  double start_angle;
  double end_angle;
  Nearest_params temp_nearest_params(Mat<3> p, Constraint is_feasible, double max_distance) const override;
  Mat<3> temp_point(Mat<2> params) const override;
  private:
  Tree_curve _tree;
};

typedef std::vector<std::unique_ptr<Entity<1>>> Composite_curve;

class Trimmed_surface : public Entity<2> {
  public:
  Trimmed_surface() : parametric_segments(n_div) {}
  bool inside(Mat<2> params) const;
  std::unique_ptr<Entity<2>> surface;
  std::vector<std::unique_ptr<Composite_curve>> curves;
  std::vector<std::vector<Mat<2>>> parametric_segments;
  Mat<3> nearest_point(Mat<3> p, double max_distance = default_max_dist) const override;
  Nearest_params temp_nearest_params(Mat<3> p, Constraint is_feasible, double max_distance) const override;
  inline Mat<3> temp_point(Mat<2> p) const override {return surface->point(p);}
};

class Geom : public Surface_geom {
  public:
  Geom(std::string file_name);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  inline std::vector<double> intersections(Mat<> point0, Mat<> point1) override {return {};}
  inline next::Sequence<Geom_edge&> edges() override {return next::Sequence<Geom_edge&>::vector_view(_edges);}
  void visualize(std::string file_name);
  private:
  std::vector<std::unique_ptr<Trimmed_surface>> _surfaces;
  std::vector<Geom_edge> _edges;
};
#endif

}
#endif
