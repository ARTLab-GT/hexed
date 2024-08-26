#ifndef HEXED_BREP_HPP_
#define HEXED_BREP_HPP_

#include <memory>
#include <functional>
#include "math.hpp"
#include "constants.hpp"
#include "Tree_curve.hpp"
#include "Surface_geom.hpp"

namespace hexed::brep {

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
    return temp_nearest_params(trans_mat.transform.colPivHouseholderQr().solve(p/scale - trans_mat.translate),
                               is_feasible, max_distance);
  }
  Mat<3> nearest_point(Mat<3> p, double max_distance = default_max_dist) const {
    Nearest_params params = nearest_params(p, [](Mat<n_param>){return true;}, max_distance);
    HEXED_ASSERT(params.is_feasible, "no feasible nearest point found");
    return point(params.params);
  }
  Mat<3> point(Mat<n_param> params) const {
    return _convert(temp_point(params));
  }
  virtual void reparameterize(Mat<n_param, 2> bounds) {}
  Int n_div = math::pow(2, 10);
  double scale = 1.;
  constexpr static double default_max_dist = std::sqrt(huge);
  Trans_mat trans_mat;
  Mat<3> _convert(Mat<3> p) const {
    return scale*(trans_mat.translate + trans_mat.transform*p);
  }
  protected:
  virtual Nearest_params temp_nearest_params(Mat<3>, Constraint is_feasible, double max_distance) const {
    return {Mat<n_param>::Zero(), false};
  };
  virtual Mat<3> temp_point(Mat<n_param>) const = 0;
  private:
};

class Circular_arc : public Entity<1> {
  public:
  inline Circular_arc(Mat<3> c, double r, double s, double e) : center{c}, radius{r}, start_angle{s}, end_angle{e} {}
  Mat<3> center;
  double radius;
  double start_angle;
  double end_angle;
  protected:
  Nearest_params temp_nearest_params(Mat<3>, Constraint is_feasible, double max_distance) const override;
  Mat<3> temp_point(Mat<1>) const override;
};

class Line_segment : public Entity<1> {
  public:
  inline Line_segment(Mat<3, 2> endpts) : endpoints{endpts} {}
  Mat<3, 2> endpoints;
  protected:
  Nearest_params temp_nearest_params(Mat<3>, Constraint is_feasible, double max_distance) const override;
  Mat<3> temp_point(Mat<1> params) const override {return endpoints*Mat<2>{1. - params(0), params(0)};}
};

class Plane : public Entity<2> {
  public:
  inline Plane(Mat<3> origin, Mat<3, 2> coord_vectors) : _origin{origin}, _vecs{coord_vectors} {}
  protected:
  void reparameterize(Mat<2, 2> bounds) override {
    _origin = temp_point(bounds(all, 0));
    _vecs = _vecs*(bounds(all, 1) - bounds(all, 0)).asDiagonal();
  }
  Nearest_params temp_nearest_params(Mat<3> p, Constraint is_feasible, double max_distance) const override {
    Mat<2> params = _vecs.colPivHouseholderQr().solve(p - _origin);
    return {params, is_feasible(params)};
  }
  inline Mat<3> temp_point(Mat<2> params) const override {return _origin + _vecs*params;}
  private:
  Mat<3> _origin;
  Mat<3, 2> _vecs;
};

class Revolution_surface : public Entity<2> {
  public:
  Revolution_surface(Entity<1>*, Line_segment, double start_angle = 0, double end_angle = 2*constants::pi);
  Mat<3> rotate(Mat<3>, double angle) const;
  std::unique_ptr<Entity<1>> generatrix;
  Line_segment axis;
  double start_angle;
  double end_angle;
  protected:
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
  protected:
  Nearest_params temp_nearest_params(Mat<3> p, Constraint is_feasible, double max_distance) const override;
  inline Mat<3> temp_point(Mat<2> p) const override {return surface->point(p);}
};

class Geom : public Surface_geom {
  public:
  Geom(std::string file_name);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = Entity<2>::default_max_dist, double distance_guess = Entity<2>::default_max_dist) override;
  inline std::vector<double> intersections(Mat<> point0, Mat<> point1) override {return {};}
  inline next::Sequence<Geom_edge&> edges() override {return next::Sequence<Geom_edge&>::vector_view(_edges);}
  void visualize(std::string file_name) const;
  private:
  std::vector<std::unique_ptr<Trimmed_surface>> _surfaces;
  std::vector<Geom_edge> _edges;
};

}
#endif
