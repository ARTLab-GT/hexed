#ifndef HEXED_SURFACE_GEOM
#define HEXED_SURFACE_GEOM

#include "math.hpp"

namespace hexed
{

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
class Surface_geom
{
  public:
  virtual ~Surface_geom() = default;
  //! computes the point on the surface which is nearest to `point`
  virtual Mat<> nearest_point(Mat<> point) = 0;
  /*! \brief Computes the set of intersection points between a line and the surface.
   * \details The line is defined parametrically to be the set of points
   * \f$ [\text{point0}] + t [\text{point1}] \f$ for all \f$ t \in \mathbb{R} \f$.
   * Returns the (potentially empty) set of \f$ t \f$ values where the line intersects the surface.
   */
  virtual std::vector<double> intersections(Mat<> point0, Mat<> point1) = 0;
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
  Mat<> nearest_point(Mat<> point) override; //!< \note `point` must not be `center`
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override;
};

}
#endif
