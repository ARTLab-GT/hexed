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

class Coordinate_change {
  public:
  Coordinate_change(Mat<3> translate = Mat<3>::Zero(), Mat<3, 3> transform = Mat<3, 3>::Identity());
  inline Mat<3> to_model(Mat<3> definition) const {return _translate + _transform*definition;}
  inline Mat<3> to_definition(Mat<3> model) const {return _inv*(model - _translate);}
  inline Mat<3, 3> transform() const {return _transform;}
  inline Mat<3> translate() const {return _translate;}
  Coordinate_change operator()(Coordinate_change that) const;
  private:
  Mat<3> _translate;
  Mat<3, 3> _transform;
  Mat<3, 3> _inv;
};

template <int n_param>
class Parametric {
  public:
  struct Nearest_parameters {
    Mat<n_param> params;
    bool is_feasible;
  };
  typedef std::function<bool(Mat<n_param>)> Constraint;
  virtual ~Parametric() = default;
  virtual Mat<3> point(Mat<n_param> params) const = 0;
  virtual Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const = 0;
  Mat<3> nearest_point(Mat<3> p) const {
    return point(nearest_params(p, [](Mat<n_param>){return true;}, default_max_dist).params);
  }
  virtual void reparameterize(Mat<n_param, 2> bounds) {}
};

template <int n_param>
class Transformed : public Parametric<n_param> {
  public:
  Transformed(Parametric<n_param>* param, Coordinate_change coord) : _param{param}, _coord{coord} {}
  Mat<3> point(Mat<n_param> params) const override {return _coord.to_model(_param->point(params));}
  Parametric<n_param>::Nearest_parameters nearest_params(
    Mat<3> p, Parametric<n_param>::Constraint is_feasible, double max_distance
  ) const override {
    return _param->nearest_params(_coord.to_definition(p), is_feasible, max_distance);
  }
  private:
  std::unique_ptr<Parametric<n_param>> _param;
  Coordinate_change _coord;
};

class Line_segment : public Parametric<1> {
  public:
  inline Line_segment(Mat<3, 2> endpoints) : _endpoints{endpoints} {}
  inline Mat<3> point(Mat<1> params) const override {return _endpoints*Mat<2>{1. - params(0), params(0)};}
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  private:
  Mat<3, 2> _endpoints;
};

class Circular_arc : public Parametric<1> {
  public:
  Circular_arc(Mat<3> center, double radius, double start_angle, double end_angle);
  Mat<3> point(Mat<1> params) const override;
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
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

class Revolution_surface : public Parametric<2> {
  public:
  Revolution_surface(Parametric<1>* generatrix, Line_segment, Int n_div, double start_angle = 0, double end_angle = 2*constants::pi);
  Mat<3> rotate(Mat<3>, double angle) const;
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  inline Mat<3> point(Mat<2> params) const override;
  private:
  class _Find_nearest;
  std::unique_ptr<Parametric<1>> _generatrix;
  Line_segment _axis;
  Int _n_div;
  double _start_angle;
  double _end_angle;
  Tree_curve _tree;
};

typedef std::vector<std::unique_ptr<Parametric<1>>> Composite_curve;

class Trimmed_surface {
  public:
  Trimmed_surface(Parametric<2>* surface, std::vector<std::unique_ptr<Composite_curve>>&& curves);
  inline const Parametric<2>& surface() const {return *_surf;}
  next::Sequence<const Composite_curve&> curves() const;
  bool inside(Mat<2> parameters) const;
  private:
  std::unique_ptr<Parametric<2>> _surf;
  std::vector<std::unique_ptr<Composite_curve>> _curves;
};

class Geom_3d {
  public:
  Geom_3d(std::string file_name, Int n_div);
  void visualize(std::string format, std::string file_name,
                 Int n_div = 100, bool vis_volume = true, Mat<3, 2> bounds = Mat<3>::Ones()*Mat<2>::Unit(1).transpose()) const;
  private:
  std::vector<Trimmed_surface> _surfaces;
};

#if 0

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
