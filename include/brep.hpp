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

  //! \brief Represents an intersection between a line (defined by 2 points) and some parametric geometry.
  struct Intersection_parameters {
    Mat<n_param> params; //!< \brief the parameters of the intersection point on the parametric geometry
    double interp_coef; //!< \brief coefficient of the intersection point as an interpolation between the 2 line points
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

  //! \brief returns the set of intersections between the infinite line passing through `points` and `this`
  //! \attention for `n_param == 1`, the calculation shall be treated as 2D and the last coordinate shall be ignored
  virtual std::vector<Intersection_parameters> intersection_params(Mat<3, 2> points) const = 0;

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

  /*! \brief whether it is necessary to check the boundary curves of a `Trimmed_surface` of this surface
   * even if a feasible nearest point was found.
   * \details The default implementation returns `true`,
   * which should be correct for most derived classes.
   */
  virtual bool must_check_boundary() const {return true;}

  //! \brief returns the parameter bounds in the original IGES definition (low, high)
  virtual Mat<2, n_param> orig_param_bounds() const = 0;
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
  typedef typename Parametric<n_param>::Intersection_parameters Inter_par;
  std::vector<Inter_par> intersection_params(Mat<3, 2> points) const override {
    for (int col = 0; col < 2; ++col) points(all, col) = _coord.to_definition(points(all, col));
    return _param->intersection_params(points);
  }
  //! \brief forwards to transformed entity
  bool must_check_boundary() const override {return _param->must_check_boundary();}
  //! \brief forwards to transformed entity
  Mat<2, n_param> orig_param_bounds() const override {return _param->orig_param_bounds();}
  private:
  std::unique_ptr<Parametric<n_param>> _param;
  Coordinate_change _coord;
};

//! \brief Line segment specified by its endpoints.
class Line_segment : public Parametric<1> {
  public:
  //! \details Each column of `endpoints` is an endpoint of the segment.
  //! They can be retrieved by `point({0.})` and `point({1.})`, respectively.
  Line_segment(Mat<3, 2> endpoints);
  inline Mat<3> point(Mat<1> params) const override {return _endpoints*Mat<2>{1. - params(0), params(0)};}
  inline double length() const {return _length;}
  //! \details Endpoints are included in nearest point search.
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  std::vector<Intersection_parameters> intersection_params(Mat<3, 2> points) const override;
  Mat<2, 1> orig_param_bounds() const override;
  private:
  Mat<3, 2> _endpoints;
  double _length;
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
  std::vector<Intersection_parameters> intersection_params(Mat<3, 2> points) const override;
  Mat<2, 1> orig_param_bounds() const override;
  private:
  Mat<3> _center;
  double _radius;
  double _start_angle;
  double _end_angle;
};

//! \brief a plane represented parametrically by an origin and coordinate vectors
class Plane : public Parametric<2> {
  public:
  //! \brief Specify `origin` and `coord_vectors` such that `point(p)` will yield `origin + coord_vectors*p`.
  inline Plane(Mat<3> origin, Mat<3, 2> coord_vectors) : _origin{origin}, _vecs{coord_vectors} {}
  //! \note Does not include boundary points in search.
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible,
                                    double max_distance) const override;
  std::vector<Intersection_parameters> intersection_params(Mat<3, 2> points) const override;
  inline Mat<3> point(Mat<2> params) const override {return _origin + _vecs*params;}
  /*! \brief Reparameterizes the plane
   * to contain the points `point(p)` where `bounds(i, 0) <= p(i) && p(i) <= bounds(i, 1)`.
   * \details `bounds(i, j)` need not be in [0, 1], but must statisfy `bounds(i, 0) < bounds(i, 1)`.
   * \see `Parametric::reparameterize`
   */
  Mat<2, 2> reparameterize(Mat<2, 2> bounds) override;
  //! \brief returns `false`; if the nearest point on the plane is feasible, there is no need to check the boundary
  inline bool must_check_boundary() const override {return false;}
  Mat<2, 2> orig_param_bounds() const override;
  private:
  Mat<3> _origin;
  Mat<3, 2> _vecs;
};

//! \brief surface constructed by revolving a 3D curve about an arbitrary axis
class Revolution_surface : public Parametric<2> {
  public:
  /*!
   * \param generatrix Curve to be revolved. Acquires ownership of `generatrix`. Need not be coplanar with `axis`.
   * \param axis Axis about which to revolve `generatrix`.
   *             The length does not matter, but the direction does,
   *             because it determines the direction of rotation by the right hand rule.
   *             Thus swapping the endpoints reverses the sense of rotation.
   * \param n_div For the purpose of nearest point calculations,
   *              `generatrix` shall be discretized into a polygonal curve with `n_div` segments,
   *              with uniform spacing in parameter space.
   *              Must be a power of 2.
   * \param start_angle Same behavior and requirements as for `Circular_arc::Circular_arc()`
   * \param end_angle Same behavior and requirements as for `Circular_arc::Circular_arc()`
   */
  Revolution_surface(Parametric<1>* generatrix, Line_segment axis, Int n_div,
                     double start_angle = 0, double end_angle = 2*constants::pi);
  //! \brief Rotates `p` about the axis by `angle`.
  //! \details Helper function made available to you cause why not?
  Mat<3> rotate(Mat<3> p, double angle) const;
  /*! \details `param(0)` specifies the point on the generatrix.
   * param(1) specifies the angle of rotation such that `param(1) = 0.` yields `start_angle`
   * and `param(1) = 1.` yields `end_angle`.
   */
  Mat<3> point(Mat<2> params) const override;
  //! \note Includes boundary points in search.
  Nearest_parameters nearest_params(Mat<3> point, Constraint is_feasible, double max_distance) const override;
  std::vector<Intersection_parameters> intersection_params(Mat<3, 2> points) const override;
  Mat<2, 2> orig_param_bounds() const override;
  private:
  class _Find_nearest;
  class _Find_intersects;
  double _unlimited_best_angle(Mat<3> arc_point, Mat<3> radius) const;
  double _limited_best_angle(Mat<3> arc_point, Mat<3> radius) const;
  Mat<3> _best_point(Mat<3> arc_point, Mat<3> radius) const;
  std::unique_ptr<Parametric<1>> _generatrix;
  Line_segment _axis;
  Mat<3> _unit_axis;
  Int _n_div;
  double _start_angle;
  double _end_angle;
  Tree_curve _tree;
};

template <int n_param>
class Nurbs : public Parametric<n_param> {
  public:
  Nurbs(std::vector<Array<double>> knots, Array<double> weights, Array<double> control_points, Int n_div,
        Mat<2, n_param> param_bounds);
  Mat<3> point(Mat<n_param> params) const override;
  Parametric<n_param>::Nearest_parameters
    nearest_params(Mat<3> point, Parametric<n_param>::Constraint is_feasible, double max_distance) const override;
  std::vector<typename Parametric<n_param>::Intersection_parameters>
    intersection_params(Mat<3, 2> points) const override;
  Mat<2, n_param> orig_param_bounds() const override;
  private:
  Int _find_knot(int i_dim, double param) const;
  struct _Nearest_params {
    Mat<3> target;
    Mat<n_param> params;
    bool is_feasible;
    double dist_sq;
    Int n_eval;
  };
  void _recursive_nearest(_Nearest_params&, std::array<Int, n_param> start_node, Int size,
                          Parametric<n_param>::Constraint is_feasible) const;
  double _max_deriv;
  std::vector<Array<double>> _knots;
  // finds the knot at the start of the interval bracketing `param` along the `i_dim`th parameter axis
  std::array<Int, n_param> _n_basis;
  std::array<int, n_param> _degree;
  Array<double> _weights;
  Array<double> _control_points;
  Int _n_div;
  Mat<2, n_param> _orig_bounds;
};

//! \brief A list of curves, where the end point of each should coincide with start of the next.
typedef std::vector<std::unique_ptr<Parametric<1>>> Composite_curve;

/*! \brief A surface created by trimming a parametric surface with closed curves.
 * \details Specifically, given a parametric surface and a set of closed curves on that surface,
 * the resulting trimmed surface is the set of points on the surface such a ray originating from that point
 * in parameter space intersects the set of closed curves an odd number of times.
 * In other words, the set of points which are inside the curves in parameter space.
 * Usually, the bounding curves are non-intersecting
 * and include one outer boundary in addition to zero or more inner boundaries,
 * but this implementation does not require that.
 * There are also no orientation requirements.
 */
class Trimmed_surface {
  public:
  /*!
   * \param surface Parametric surface to be trimmed. Acquires ownership of `surface`.
   * \param curves Bounding curves _in model space_, not parameter space.
   *               These curves should (approximately) lie on the surface.
   *               Any deviation from the surface will be a source of numerical error.
   *               \todo Implement another constructor that accepts curves in parameter space.
   * \param n_div For some calculations, the bounding curves will be discretized
   *              into polygonal curves with O(`n_div`) segments.
   *              Must be a power of 2.
   */
  Trimmed_surface(Parametric<2>* surface, std::vector<Composite_curve>&& curves,
                  std::vector<bool> is_model_space, Int n_div);
  //! \brief Access the parametric surface.
  inline const Parametric<2>& surface() const {return *_surf;}
  //! \brief Test whether a point `parameters` is inside the bounding curves in parameter space.
  bool is_inside(Mat<2> parameters) const;
  //! \brief Access the bounding curves.
  next::Sequence<const Tree_curve&> curves() const;
  /*! \brief Compute the point on the trimmed surface (including the boundary) nearest to `point`.
   * \details If the nearest point would be further than `max_dist` from `point`,
   * the empty `Nearest_point` is returned.
   */
  Nearest_point<3> nearest_point(Mat<3> point, double max_dist) const;
  std::vector<double> intersections(Mat<3, 2> endpoints) const;
  private:
  // Performs the real initialization work once the curves have been discretized.
  // Discretization is performed by the constructor.
  void _initialize(std::vector<std::vector<std::vector<Mat<2>>>>& curves);
  Int _n_div;
  double _sz;
  std::unique_ptr<Parametric<2>> _surf;
  std::vector<Tree_curve> _curves;
  // No simple way to explain this.
  // Need to write a dedicated article about distinguishing inside/outside points, which _param_segmetns is a part of.
  std::vector<std::vector<Mat<2>>> _param_segments;
};

//! \brief a `Surface_geom` consisting of a set of `Parametric<1>` curves
class Geom_2d : public Surface_geom {
  public:
  //! \param file_name Name of file containing geometry. Must be in IGES format.
  //! \param n_div Any entities that need to be discretized will be so with `n_div` subdivisions. Must be a power of 2.
  Geom_2d(std::string file_name, Int n_div);
  /*! \brief Writes visualization files of the geometry to help diagnose import/translation bugs.
   * \details For visualization purposes, entities will be discretized with `n_div` segments.
   * This is not the same as the `n_div` passed to the constructor, and need not be a power of 2.
   * It should usually be much less than the `n_div` passed to the constructor,
   * because visualization is more expensive than nearest-point calculations and requires less precision.
   * If the input file is named `INPUT_FILE`
   * and the file extension of the specified visualization format is `EXT`, this function
   * will write a file `INPUT_FILE_curves.EXT` with all curves in the geometry.
   * Even though `Geom_2d` is supposed to represent a 2D geometry, the curves will be 3D.
   * This is because the underlying representation is 3D,
   * and if this turns out not to lie in the \f$ x_0, x_1 \f$ plane,
   * this is a potential source of problems and important information to convey in the visualization file.
   */
  void visualize(std::string format, std::string file_name, Int n_div = 100);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  //! \brief Dummy implementation that returns an empty vector.
  inline std::vector<double> intersections(Mat<> point0, Mat<> point1) override {return {};}
  next::Sequence<Mat<3>> points() override;
  private:
  std::vector<std::unique_ptr<Parametric<1>>> _curves;
};

//! \brief a `Surface_geom` consisting of a set of 3D `Trimmed_surface`s
class Geom_3d : public Surface_geom {
  public:
  //! \param file_name Name of file containing geometry. Must be in IGES format.
  //! \param n_div Any entities that need to be discretized will be so with `n_div` subdivisions. Must be a power of 2.
  Geom_3d(std::string file_name, Int n_div);
  /*! \brief Writes visualization files of the geometry to help diagnose import/translation bugs.
   * \details For visualization purposes, entities will be discretized with `n_div` segments.
   * This is not the same as the `n_div` passed to the constructor, and need not be a power of 2.
   * It should usually be much less than the `n_div` passed to the constructor,
   * because visualization is more expensive than nearest-point calculations and requires less precision.
   *
   * If the input file is named `INPUT_FILE`
   * and the file extension of the specified visualization format is `EXT`, the visualization files are:
   * - `FILE_NAME_curves.EXT`: Contains all bounding curves of all surfaces.
   * - `FILE_NAME_surfaces.EXT`: Contains all parametric surfaces
   *   (for all parameters in \f$ [0, 1] \times [0, 1] \f$).
   *   The surfaces have a field variable "inside"
   *   which is set to 1 for all points that are inside the bounding curves and 0 for all that are outside.
   *   Pro tip: In Paraview, you can get a sense of the actual geometry by enabling opacity mapping
   *   so that the trimmed regions are transparent.
   * - `FILE_NAME_distance.EXT` (only if `vis_volume = true`): A 3D block with corners given by the columns of `bounds`
   *   and a field variable indicating the distance from the nearest point on the surface.
   *   This can be a good way to debug the distance calculations, assuming the geometry is correct.
   */
  void visualize(std::string format, std::string file_name,
                 Int n_div = 100, bool vis_volume = true, Mat<3, 2> bounds = Mat<3>::Ones()*Mat<2>::Unit(1).transpose());
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
  next::Sequence<const Tree_curve&> edges() override;
  private:
  std::vector<Trimmed_surface> _surfaces;
};

}
#endif
