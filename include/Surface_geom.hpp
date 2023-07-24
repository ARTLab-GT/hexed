#ifndef HEXED_SURFACE_GEOM
#define HEXED_SURFACE_GEOM

#include <memory>
#include "math.hpp"

namespace hexed
{

/*! \brief Represents a surface geometry implicitly for meshing.
 * \details Abstract class which represents geometry by supporting
 * the operations `Surface_geom::nearest_point`/`Surface_geom::bounded_nearest` and `Surface_geom::intersections`.
 * Either of these operations fully defines the surface geometry
 * and each could, in theory, be defined in terms of the other.
 * However, in practice they are defined separately
 * so that derived classes can implement each in the most efficient and robust way possible.
 * It is the responsibility of the derived class to ensure
 * that they are mutually consistent.
 * The derived class may impose dimensionality restrictions on the arguments.
 * Implementations must be thread-safe.
 * \note Depending on whether `nearest_point` or `bounded_nearest` is easier to implement efficiently,
 * consider subclassing `Bounded_search_geom` or `Unbounded_search_geom`.
 */
class Surface_geom
{
  public:
  virtual ~Surface_geom() = default;
  /*! \brief Computes the point on the surface which is nearest to `point`.
   * \param distance_guess If you have some reason to suspect the nearest point is within a certain distance
   * of the input point, you can pass it to this parameter as a hint to improve performance.
   */
  virtual Mat<> nearest_point(Mat<> point, double distance_guess = std::numeric_limits<double>::max()) = 0;
  /*! \brief Computes the nearest point, if any, within a certain radius.
   * \details If there is a point on the geometry within `max_distance` of `point`,
   * you will get it as a `math::Nearest_point` object.
   * If there is no such point, you will get an empty `math::Nearest_point`.
   */
  virtual math::Nearest_point<dyn> bounded_nearest(Mat<> point, double max_distance) = 0;
  /*! \brief Computes the set of intersection points between a line and the surface.
   * \details The line is defined parametrically to be the set of points
   * \f$ [\text{point0}] + t [\text{point1}] \f$ for all \f$ t \in \mathbb{R} \f$.
   * Returns the (potentially empty) set of \f$ t \f$ values where the line intersects the surface.
   */
  virtual std::vector<double> intersections(Mat<> point0, Mat<> point1) = 0;
};

/*! \brief A `Surface_geom` where `nearest_point` is implemented in terms of `bounded_nearest`.
 * \details If you want to make a `Surface_geom` subclass
 * where searching within a finite radius is inherently easier to implement efficiently,
 * this class can provide a convenient and reasonably fast generalization to unbounded searches.
 */
class Bounded_search_geom
{
  public:
  /*! \details Calls `bounded_nearest(point, distance_guess)` and if no point is found,
   * iteratively doubles the search radius until it finds a point.
   * If the speed of your `bounded_nearest` implementation is strongly dependent on the search radius,
   * this should be a pretty fast way to compute `nearest_point`.
   */
  Mat<> nearest_point(Mat<> point, double distance_guess = std::numeric_limits<double>::max()) override;
};

/*! \brief A `Surface_geom` where `bounded_nearest` is implemented in terms of `nearest_point`.
 * \details If you want to make a `Surface_geom` subclass
 * where there is no particular benefit to searching within a finite radius,
 * this class can provide a convenient generalization to bounded searches.
 */
class Unbounded_search_geom
{
  public:
  //! \details calls `nearest_point(point, max_distance)` and returns an empty point if it is too far away.
  math::Nearest_point<dyn> bounded_nearest(Mat<> point, double max_distance) override;
};

/*! \brief Combines multiple `Surface_geom`s into one.
 * \details The nearest point to the compound geometry
 * is the nearest of the nearest points on the component geometries
 * and the intersection set is the union of the intersection sets of the components.
 */
class Compound_geom : public Surface_geom
{
  std::vector<std::unique_ptr<Surface_geom>> components;
  public:
  Compound_geom(std::vector<Surface_geom*>); //!< acquires ownership
  math::Nearest_point<dyn> bounded_nearest(Mat<> point, double max_distance) override;
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
};

/*! \brief Represents hypersphere in any dimensionality.
 * \details Can behave as an interval, cylinder, sphere, or higher-dimensional equivalent
 * depending on the number of entries in the point vectors supplied as arguments.
 * However, all points supplied must have the same dimensionality,
 * or else behavior is undefined.
 */
class Hypersphere : public Surface_geom
{
  Mat<> c;
  double r;
  public:
  Hypersphere(Mat<> center, double radius);
  math::Nearest_point<dyn> bounded_nearest(Mat<> point, double max_distance) override; //!< \note `point` must not be `center`
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
};

}
#endif
