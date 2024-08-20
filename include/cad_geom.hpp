#ifndef HEXED_CAD_GEOM_HPP_
#define HEXED_CAD_GEOM_HPP_

#include <memory>
#include "math.hpp"
#include "constants.hpp"

namespace hexed::cad_geom {

struct Trans_mat {
  Mat<3, 3> transform = Mat<3, 3>::Identity();
  Mat<3> translate = Mat<3>::Zero();
};

template <int n_param>
class Entity {
  public:
  virtual ~Entity() = default;
  virtual Mat<3> temp_nearest_point(Mat<3>) const {return Mat<3>::Zero();};
  virtual Mat<3> temp_point(Mat<n_param>) const = 0;
  Mat<3> nearest_point(Mat<3> p) const {
    return _convert(temp_nearest_point(trans_mat.transform.colPivHouseholderQr().solve(p/scale - trans_mat.translate)));
  }
  Mat<3> point(Mat<n_param> params) const {
    return _convert(temp_point(params));
  }
  int n_div = math::pow(10, 2);
  double scale = 1.;
  Trans_mat trans_mat;
  private:
  Mat<3> _convert(Mat<3> p) const {
    return scale*(trans_mat.translate + trans_mat.transform*p);
  }
};

class Circular_arc : public Entity<1> {
  public:
  inline Circular_arc(Mat<3> c, double r, double s, double e) : center{c}, radius{r}, start_angle{s}, end_angle{e} {}
  Mat<3> temp_point(Mat<1>) const override;
  Mat<3> center;
  double radius;
  double start_angle;
  double end_angle;
};

class Line_segment : public Entity<1> {
  public:
  inline Line_segment(Mat<3, 2> endpts) : endpoints{endpts} {}
  Mat<3> temp_point(Mat<1> params) const override {return endpoints*Mat<2>{1. - params(0), params(0)};}
  Mat<3, 2> endpoints;
};

class Plane : public Entity<2> {
  public:
  inline Plane(Mat<3> origin, Mat<3, 2> coord_vectors) : _origin{origin}, _vecs{coord_vectors} {}
  Mat<3> temp_nearest_point(Mat<3> p) const override {
    return _origin + _vecs*_vecs.colPivHouseholderQr().solve(p - _origin);
  }
  inline Mat<3> temp_point(Mat<2> params) const override {return _origin + _vecs*params;}
  private:
  Mat<3> _origin;
  Mat<3, 2> _vecs;
};

class Revolution_surface : public Entity<2> {
  public:
  Revolution_surface(Entity<1>*, Line_segment, double start_angle = 0, double end_angle = 2*constants::pi);
  Mat<3> temp_nearest_point(Mat<3> p) const override;
  Mat<3> temp_point(Mat<2> params) const override;
  std::unique_ptr<Entity<1>> generatrix;
  Line_segment axis;
  double start_angle;
  double end_angle;
};

typedef std::vector<std::unique_ptr<Entity<1>>> Composite_curve;

class Trimmed_surface : public Entity<2> {
  public:
  inline Mat<3> temp_nearest_point(Mat<3> p) const override{return surface->nearest_point(p);};
  inline Mat<3> temp_point(Mat<2> p) const override {return surface->point(p);}
  std::unique_ptr<Entity<2>> surface;
  std::vector<std::unique_ptr<Composite_curve>> curves;
};

class Geom {
  public:
  Geom(std::string file_name);
  void visualize(std::string file_name) const;

  private:
  std::vector<std::unique_ptr<Trimmed_surface>> _surfaces;
};

}
#endif
