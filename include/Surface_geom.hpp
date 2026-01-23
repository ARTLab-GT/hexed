#ifndef HEXED_SURFACE_GEOM_HPP_
#define HEXED_SURFACE_GEOM_HPP_

#include <memory>
#include "math.hpp"
#include "Nearest_point.hpp"
#include "Sequence.hpp"
#include "Tree_curve.hpp"

namespace hexed {

class Geom_edge {
  public:
  virtual Mat<3> point(double) const = 0;
  virtual double arc_length(double) const = 0;
  virtual Mat<3> tangent_average(double) const = 0;
  virtual double tangent_radius(double) const = 0;
  virtual double arg_nearest_point(Mat<3>) const = 0;
};

class Compound_edge : public Geom_edge {
  public:
  Compound_edge(std::vector<std::shared_ptr<Geom_edge>> _edges, std::vector<bool> reverse);
  Mat<3> point(double) const override;
  double arc_length(double) const override;
  Mat<3> tangent_average(double) const override;
  double tangent_radius(double) const override;
  double arg_nearest_point(Mat<3>) const override;
  private:
  std::pair<Int, double> _global_index(double) const;
  std::vector<std::shared_ptr<Geom_edge>> _edges;
  std::vector<bool> _reverse;
  std::vector<double> _arc_length_start;
};

/*! \brief Represents a surface geometry implicitly for meshing.
 * \details Abstract class which represents geometry by supporting
 * the operations `Surface_geom::nearest_point` and `Surface_geom::intersections`.
 * Either of these operations fully defines the surface geometry
 * and each could, in theory, be defined in terms of the other.
 * However, in practice they are defined separately
 * so that derived classes can implement each in the most efficient and robust way possible.
 * It is the responsibility of the derived class to ensure
 * that they are mutually consistent.
 * The derived class may impose dimensionality restrictions on the arguments.
 * Implementations must be thread-safe.
 */
class Surface_geom {
  public:
  virtual ~Surface_geom() = default;
  /*! \brief Computes the point on the surface which is nearest to `point` within `max_distance`.
   * \details If no point is found, returns an empty `Nearest_point`.
   * \param point The point you want to find the nearest point to.
   * \param max_distance Only consider points within `max_distance` of `point`.
   * \param distance_guess If you have some reason to suspect the nearest point is within a certain distance
   * of the input point, you can pass it to this parameter as a hint to possibly improve performance.
   */
  virtual Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) = 0;
  /*! \brief Computes the set of intersection points between a line and the surface.
   * \details The line is defined parametrically to be the set of points
   * \f$ [\text{point0}] + t [\text{point1}] \f$ for all \f$ t \in \mathbb{R} \f$.
   * Returns the (potentially empty) set of \f$ t \f$ values where the line intersects the surface.
   * \deprecated Newer methods should be based purely on nearest point projections.
   */
  virtual std::vector<double> intersections(Mat<> point0, Mat<> point1, bool high_precision = true) = 0;
  /*! \brief Returns a list of any geometry edges that require mesh edges to be snapped to them.
   * \details Only used in 3D.
   * Default implementation returns an empty sequence, but derived classes may override.
   */
  inline virtual next::Sequence<const Geom_edge&> edges() {return {};}
  inline virtual next::Sequence<Mat<3>> points() {return {};}
};

/*! \brief Combines multiple `Surface_geom`s into one.
 * \details The nearest point to the compound geometry
 * is the nearest of the nearest points on the component geometries
 * and the intersection set is the union of the intersection sets of the components.
 */
class Compound_geom : public Surface_geom {
  public:
  Compound_geom(std::vector<Surface_geom*>); //!< acquires ownership
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  std::vector<double> intersections(Mat<> point0, Mat<> point1, bool high_precision = true) override;
  next::Sequence<const Geom_edge&> edges() override;
  next::Sequence<Mat<3>> points() override;
  private:
  std::vector<std::unique_ptr<Surface_geom>> _components;
};

/*! \brief Represents hypersphere in any dimensionality.
 * \details Can behave as an interval, cylinder, sphere, or higher-dimensional equivalent
 * depending on the number of entries in the point vectors supplied as arguments.
 * However, all points supplied must have the same dimensionality,
 * or else behavior is undefined.
 */
class Hypersphere : public Surface_geom {
  public:
  Hypersphere(Mat<> center, double radius);
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override;
  std::vector<double> intersections(Mat<> point0, Mat<> point1, bool high_precision = true) override;
  private:
  Mat<> _center;
  double _radius;
};

}
#endif
