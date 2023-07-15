#ifndef HEXED_SIMPLEX_GEOM_HPP_
#define HEXED_SIMPLEX_GEOM_HPP_

#include "Surface_geom.hpp"
#include "global_hacks.hpp"

namespace hexed
{

/*! \brief Represents discrete geometry composed of [simplices](https://en.wikipedia.org/wiki/Simplex).
 * \details This can be used as an interface for geometry derived from an STL file in 3D
 * or a list of node coordinates in 2D.
 * Each simplex is represented as an `n_dim` by `n_dim` matrix where each column is the coordinates of one vertex.
 * The order of the vertices is arbitrary, and there are no requirements on inter-simplex continuity.
 * Since the `Surface_geom` interface is so minimal, the geometry is simply viewed as a collection of unrelated simplices,
 * so we do not care about orientation or watertightness.
 * The you are free to modify the simplex list at will, since there are no requirements on it.
 * All input points must have exactly `n_dim` entries.
 */
template <int n_dim>
class Simplex_geom : public Surface_geom
{
  public:
  std::vector<Mat<n_dim, n_dim>> simplices;
  Simplex_geom(const std::vector<Mat<n_dim, n_dim>>&  sims) : simplices{sims} {}
  Simplex_geom(      std::vector<Mat<n_dim, n_dim>>&& sims) : simplices{sims} {}

  /*! \details Iterates through all simplices and finds the nearest point on each,
   * whether that point lies in the interior or on the edge or a vertex.
   * Returns the global nearest of all those points.
   */
  Mat<> nearest_point(Mat<> point) override;

  /*! \details Evaluates intersections with each individual element and assembles a global list.
   * Intersections with the boundary of a simplex are always considered valid intersections,
   * so if the line passes exactly through the shared boundary of multiple simplices
   * then duplicate intersections may be obtained.
   */
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override
  {
    Stopwatch sw;
    sw.start();
    std::vector<double> inters;
    Mat<n_dim> diff = point1 - point0;
    for (Mat<n_dim, n_dim> sim : simplices) {
      // compute intersection between line and plane of simplex
      Mat<n_dim, n_dim> lhs;
      lhs(all, 0) = -diff;
      for (int col = 1; col < n_dim; ++col) lhs(all, col) = sim(all, col) - sim(all, 0);
      Mat<n_dim> soln = lhs.lu().solve(point0 - sim(all, 0));
      // if intersection is inside simplex, add it to the list
      Eigen::Array<double, n_dim - 1, 1> arr = soln(Eigen::seqN(1, n_dim - 1)).array();
      if ((arr >= 0.).all() && arr.sum() <= 1.) inters.push_back(soln(0));
    }
    sw.pause();
    #pragma omp atomic update
    global_hacks::numbers[1] += sw.time();
    return inters;
  }
};

typedef Simplex_geom<2> Simplex_geom2; //!< typedef for cython bindings
typedef Simplex_geom<3> Simplex_geom3; //!< typedef for cython bindings

template<> Mat<> Simplex_geom<2>::nearest_point(Mat<> point);
template<> Mat<> Simplex_geom<3>::nearest_point(Mat<> point);

/*! \brief Creates simplices from an ordered array of points representing a polygonal curve.
 * \details `points` must have exactly 2 rows.
 * Each column of `points` is the coordinates of one point.
 * The curve may or may not be closed.
 * Each point represents a unique vertex -- duplicate points are allowed,
 * but you don't need to include each interior point twice.
 * The result can be used to construct a `Simplex_geom<2>`.
 * \relates Simplex_geom
 */
std::vector<Mat<2, 2>> segments(const Mat<dyn, dyn>& points);

}
#endif
