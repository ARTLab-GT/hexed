#ifndef HEXED_BREP_HPP_
#define HEXED_BREP_HPP_

#include <memory>
#include <functional>
#include "math.hpp"
#include "constants.hpp"
#include "Tree_curve.hpp"
#include "Surface_geom.hpp"

/*! \brief Implementation of the %Boundary Representation (BRep) for solid geometry.
 * \details This namespace provides the facilities for reading BRep geometry from a file
 * and translating it into a `Surface_geom` (2D or 3D).
 * It does not provide general facilities for constructing and manipulating BRep geometry
 * (you'll need a different software package for that).
 * The main classes of interest here are `Geom_2d` and `Geom_3d`.
 * The others are helper classes used in the implementation of the `Geom_2d` and `Geom_3d`.
 */
namespace hexed::brep {

//! \brief Default maximum distance to use in nearest point calculations.
//! \details This is the largest value that can be squared while remaining finite.
const double default_max_dist = std::sqrt(huge);

//! \brief A change of coordinates from "definition" space (domain) to "model" space (codomain).
//! \details Consists of a linear transformation (usually unitary) followed by a translation.
class Coordinate_change {
  public:
  //! \note `transform` must be invertible.
  Coordinate_change(Mat<3> translate = Mat<3>::Zero(), Mat<3, 3> transform = Mat<3, 3>::Identity());
  //! \brief Applies the coordinate transformation.
  //! \details Specifically, the result is `translate() + transform()*defintion`.
  inline Mat<3> to_model(Mat<3> definition) const {return _translate + _transform*definition;}
  //! \brief Inverts the coordinate transformation.
  //! \details `to_definition(to_model(x)) == to_model(to_definition(x)) = x`.
  inline Mat<3> to_definition(Mat<3> model) const {return _inv*(model - _translate);}
  //! \brief Obtains the transformation matrix, if you care.
  inline Mat<3, 3> transform() const {return _transform;}
  //! \brief Obtains the translation vector, if you care.
  inline Mat<3> translate() const {return _translate;}
  //! \brief Returns a `Coordinate_change` which is the composition of `this` and `that`.
  //! \details I.e., `(*this)(that).to_model(x) == this->to_model(that.to_model(x))`.
  Coordinate_change operator()(Coordinate_change that) const;
  private:
  Mat<3> _translate;
  Mat<3, 3> _transform;
  Mat<3, 3> _inv;
};

/*! \brief Represents a basic geometric entity in parametric form.
 * \details `n_param` is the number of parameters the entity takes.
 * Thus `Parametric<1>` is a parametric curve and `Parametric<2>` is a parametric surface.
 * For simplicity, the number of output variables is always 3.
 * If the geometry is intended to be 2D, then the last coordinate should be (approximately) 0.
 */
template <int n_param>
class Parametric {
  public:
  //! \brief Represents the result of a `nearest_params()` calculation,
  //! \details which may or may not have found a feasible result.
  struct Nearest_parameters {
    Mat<n_param> params; //!< \brief the parameters of the nearest point on the entity, if found
    bool is_feasible; //!< \brief `true` iff a feasible nearest point was found
  };

  //! \brief Represents a constraint function for a `nearst_params()` calculation.
  typedef std::function<bool(Mat<n_param>)> Constraint;

  virtual ~Parametric() = default;

  //! \brief Obtains the point at parameters `params`.
  virtual Mat<3> point(Mat<n_param> params) const = 0;

  /*! \brief Finds the parameters of the point on the entity nearest to `point`,
   * subject to a constraint, if that point is in the interior.
   * \details This is designed as a helper function to a robust and computationally efficient implementation of
   * `Trimmed_surface::nearest_point`, so its behavior might be confusing if considered in isolation.
   * Derived classes must implement this member function.
   *
   * Consider the problem of finding the parameter vector `param`
   * to minimize the distance `(point - this->point(param)).norm()`
   * subject to the constraint that `is_feasible(param)` must be `true`
   * and the distance must be less than `max_distance`.
   * If the solution to this problem is in the interior of the entity
   * (`0 < param(i) && param(i) < 1` for all `i` in [0, `n_param`)),
   * then the implementation must return it as a `Nearest_parameters`
   * with `Nearest_parameters::is_feasible` set to `true`.
   * If there is no solution, it must return a `Nearest_parameters` with arbitrary `params`
   * and `Nearest_parameters::is_feasible` set to `false`.
   * Implementations may optionally include boundary points in their search.
   */
  virtual Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const = 0;

  /*! \brief Finds the nearest point on the entity to `p`, at least if that point is on the interior.
   * \details Wrapper for `nearest_params()` with no constraints applied.
   * \warning This is mostly for testing.
   * Due to the implementor-defined behavior for boundary points, it may not always give you what you expect.
   */
  Mat<3> nearest_point(Mat<3> p) const {
    return point(nearest_params(p, [](Mat<n_param>){return true;}, default_max_dist).params);
  }

  /*! \brief May reparameterize the entity to keep parameters in [0, 1].
   * \details For infinite entities (e.g., `Plane`),
   * it is not always known _a priori_ what parameterization makes sense,
   * so they are initialized with an arbitrary parameterization.
   * Once the range of parameters for points of interest is known,
   * this function should be called with `bounds` such that `bounds(i, 0) <= param(i) && param(i) <= bounds(i, 1)`
   * for all `i` and all points `param` that will be required in calculations.
   * This function should then modify the parameterization such that `0 <= param(i) && param(i) <= 1`
   * for all points `param` in future calculations.
   * It will then return the new bounds such that a linear function (\f$ f(x) = ax + b \f$)
   * which maps `bounds` to `this->reparameterize(bounds)`
   * will map points in the old parameterization to the new parameterization.
   * In other words, it returns the new parameter bounds in the old parameterization.
   * The default implementation does nothing and returns \f$ [\vec{0}\ \vec{1}] \f$,
   * and finite entities need not override this.
   */
  virtual Mat<n_param, 2> reparameterize(Mat<n_param, 2> bounds) {
    bounds(all, 0).setZero();
    bounds(all, 1).setOnes();
    return bounds;
  }
};

/*! \brief A parametric entity obtained by applying a `Coordinate_change` to another parametric entity
 * \details For performance and type transparency,
 * whenever possible the transformation should be applied
 * to the defining parameters of the entity (e.g. line segment endpoints) instead of using this class.
 * This class provides a solution for entities where that is not possible, such as `Circular_arc`.
 */
template <int n_param>
class Transformed : public Parametric<n_param> {
  public:
  //! \note Acquires ownership of `param`!
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

//! \brief Line segment specified by its endpoints.
class Line_segment : public Parametric<1> {
  public:
  //! \details Each column of `endpoints` is an endpoint of the segment.
  //! They can be retrieved by `point({0.})` and `point({1.})`, respectively.
  inline Line_segment(Mat<3, 2> endpoints) : _endpoints{endpoints} {}
  inline Mat<3> point(Mat<1> params) const override {return _endpoints*Mat<2>{1. - params(0), params(0)};}
  //! \details Endpoints are included in nearest point search.
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  private:
  Mat<3, 2> _endpoints;
};

//! \brief A circular arc in the \f$ x_0, x_1 \f$ plane.
class Circular_arc : public Parametric<1> {
  public:
  /*! \details The arc will contain the angle interval [`start_angle`, `end_angle`].
   * `start_angle` and `end_angle` may be in any domain
   * (i.e. they don't need to be in \f$ [-\pi, \pi] \f$ or something like that),
   * but must satisfy `end_angle` > `start_angle` (note the strict inequality).
   * If `end_angle - start_angle` \f$ \ge 2\pi \f$, then the arc will contain redundant points.
   * In this case, `nearest_params()` is not guaranteed to check more than 1 of the redundant points for feasibility,
   * and which one it checks it not specified.
   */
  Circular_arc(Mat<3> center, double radius, double start_angle, double end_angle);
  Mat<3> point(Mat<1> params) const override;
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  private:
  Mat<3> _center;
  double _radius;
  double _start_angle;
  double _end_angle;
};

//! \brief A plane represented parametrically by an origin and coordinate vectors.
class Plane : public Parametric<2> {
  public:
  //! \brief Specify `origin` and `coord_vectors` such that `point(p)` will yield `origin + coord_vectors*p`.
  inline Plane(Mat<3> origin, Mat<3, 2> coord_vectors) : _origin{origin}, _vecs{coord_vectors} {}
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                    double max_distance) const override;
  inline Mat<3> point(Mat<2> params) const override {return _origin + _vecs*params;}
  /*! \brief Reparameterizes the plane
   * to contain the points `point(p)` where `bounds(i, 0) <= p(i) && p(i) <= bounds(i, 1)`.
   * \details `bounds(i, j)` need not be in [0, 1], but must statisfy `bounds(i, 0) < bounds(i, 1)`.
   * \see `Parametric::reparameterize`
   */
  Mat<2, 2> reparameterize(Mat<2, 2> bounds) override;
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
  Trimmed_surface(Parametric<2>* surface, std::vector<Composite_curve>&& curves, Int n_div);
  inline const Parametric<2>& surface() const {return *_surf;}
  bool is_inside(Mat<2> parameters) const;
  next::Sequence<const Tree_curve&> curves() const;
  Nearest_point<3> nearest_point(Mat<3> point, double max_dist) const;
  private:
  void initialize(std::vector<std::vector<Mat<2>>>& curves);
  Int _n_div;
  double _sz;
  std::unique_ptr<Parametric<2>> _surf;
  std::vector<Tree_curve> _curves;
  std::vector<std::vector<Mat<2>>> _param_segments;
};

class Geom_3d : public Surface_geom {
  public:
  Geom_3d(std::string file_name, Int n_div);
  void visualize(std::string format, std::string file_name,
                 Int n_div = 100, bool vis_volume = true, Mat<3, 2> bounds = Mat<3>::Ones()*Mat<2>::Unit(1).transpose());
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  inline std::vector<double> intersections(Mat<> point0, Mat<> point1) override {return {};}
  next::Sequence<const Tree_curve&> edges() override;
  private:
  std::vector<Trimmed_surface> _surfaces;
};

class Geom_2d : public Surface_geom {
  public:
  Geom_2d(std::string file_name, Int n_div);
  void visualize(std::string format, std::string file_name, Int n_div = 100);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  inline std::vector<double> intersections(Mat<> point0, Mat<> point1) override {return {};}
  next::Sequence<Mat<3>> points() override;
  private:
  std::vector<std::unique_ptr<Parametric<1>>> _curves;
};

}
#endif
