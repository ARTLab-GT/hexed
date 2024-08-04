#ifndef HEXED_CAD_GEOM_HPP_
#define HEXED_CAD_GEOM_HPP_

#include <memory>
#include "math.hpp"

namespace hexed::cad_geom {

struct Trans_mat {
  Mat<3, 3> transform = Mat<3, 3>::Identity();
  Mat<3> translate = Mat<3>::Zero();
};

template <int n_param>
class Entity {
  public:
  virtual ~Entity() = default;
  virtual Mat<3> temp_point(Mat<n_param>) const = 0;
  Mat<3> point(Mat<n_param> params) const {
    return scale*(trans_mat.translate + trans_mat.transform*temp_point(params));
  }
  double scale = 1.;
  Trans_mat trans_mat;
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

class Geom {
  public:
  Geom(std::string file_name);
  void visualize(std::string file_name) const;

  private:
  std::vector<std::unique_ptr<Entity<1>>> _curves;
};

}
#endif
