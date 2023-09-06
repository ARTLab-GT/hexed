#ifndef HEXED_SIMPLEX_GEOM_HPP_
#define HEXED_SIMPLEX_GEOM_HPP_

#include "Surface_geom.hpp"
#include "Stopwatch_tree.hpp"

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
 * \see `Occt::segments`
 * \see `Occt::triangles`
 */
template <int n_dim>
class Simplex_geom : public Surface_geom
{
  #if HEXED_OBSESSIVE_TIMING
  static Stopwatch_tree stopwatch; // for benchmarking projection and intersection calculation
  #endif
  void merge(Nearest_point<n_dim>& nearest, Mat<n_dim, n_dim> sim, Mat<n_dim> point); // helper for `nearest_point`

  public:
  double gap_tol = 1e-2; //!< extend triangles by this amount, relative to their original size, to fill gaps
  std::vector<Mat<n_dim, n_dim>> simplices;
  Simplex_geom(const std::vector<Mat<n_dim, n_dim>>&  sims) : simplices{sims} {}
  Simplex_geom(      std::vector<Mat<n_dim, n_dim>>&& sims) : simplices{sims} {}
  void visualize(std::string file_name); //!< writes geometry to a Tecplot file with the specified name + file extension

  /*! \details Iterates through all simplices and finds the nearest point on each,
   * whether that point lies in the interior or on the edge or a vertex.
   * Returns the global nearest of all those points.
   * `distance_guess` is taken into account, and it is generally preferable to overestimate rather than underestimate.
   */
  Nearest_point<dyn> nearest_point(Mat<> point, double max_distance = huge, double distance_guess = huge) override
  {
    HEXED_ASSERT(distance_guess > 0, "`distance_guess` must be positive");
    #if HEXED_OBSESSIVE_TIMING
    Stopwatch sw;
    sw.start();
    #endif
    // try to find a point within `distance_guess`
    double limit = std::min(max_distance, distance_guess);
    Nearest_point<n_dim> nearest(point, limit);
    for (Mat<n_dim, n_dim> sim : simplices) {
      auto bball = math::bounding_ball(sim);
      if ((bball.center - point).norm() - std::sqrt(bball.radius_sq) <= limit) { // only bother if at least the bounding ball is close enough
        merge(nearest, sim, point);
      }
    }
    #if HEXED_OBSESSIVE_TIMING
    sw.pause();
    stopwatch.stopwatch += sw;
    stopwatch.children.at("nearest_point").stopwatch += sw;
    #pragma omp atomic update
    ++stopwatch.children.at("nearest_point").work_units_completed;
    #endif
    // if the we failed, recursively double `distance_guess` until we hit `max_distance`
    return (nearest.empty() && distance_guess < max_distance) ? nearest_point(point, max_distance, 2*distance_guess) : Nearest_point<dyn>(nearest);
  }

  /*! \details Evaluates intersections with each individual element and assembles a global list.
   * Intersections with the boundary of a simplex are always considered valid intersections,
   * so if the line passes exactly through the shared boundary of multiple simplices
   * then duplicate intersections may be obtained.
   */
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override
  {
    #if HEXED_OBSESSIVE_TIMING
    Stopwatch sw;
    sw.start();
    #endif
    std::vector<double> inters;
    Mat<n_dim> diff = point1 - point0;
    for (Mat<n_dim, n_dim> sim : simplices) {
      if (math::intersects<n_dim>(math::bounding_ball(sim), point0, point1)) { // if the line doesn't even intersect the bounding ball, don't bother
        // compute intersection between line and plane of simplex
        Mat<n_dim, n_dim> lhs;
        lhs(all, 0) = -diff;
        for (int col = 1; col < n_dim; ++col) lhs(all, col) = sim(all, col) - sim(all, 0);
        Mat<n_dim> soln = lhs.lu().solve(point0 - sim(all, 0));
        // if intersection is inside simplex, add it to the list
        Eigen::Array<double, n_dim - 1, 1> arr = soln(Eigen::seqN(1, n_dim - 1)).array();
        if ((arr >= -gap_tol).all() && arr.sum() <= 1 + gap_tol) inters.push_back(soln(0));
      }
    }
    #if HEXED_OBSESSIVE_TIMING
    sw.pause();
    stopwatch.stopwatch += sw;
    stopwatch.children.at("intersections").stopwatch += sw;
    #pragma omp atomic update
    ++stopwatch.children.at("intersections").work_units_completed;
    #endif
    return inters;
  }

  //! if compiled with `HEXED_OBSESSIVE_TIMING ON`, return a performance report. Otherwise, empty string.
  static std::string performance_report()
  {
    #if HEXED_OBSESSIVE_TIMING
    return stopwatch.report();
    #else
    return "";
    #endif
  }
};

#if HEXED_OBSESSIVE_TIMING
template <int n_dim>
Stopwatch_tree Simplex_geom<n_dim>::stopwatch("", {{"nearest_point", Stopwatch_tree("projection")}, {"intersections", Stopwatch_tree("intersection")}}); // for benchmarking projection and intersection calculation
#endif

template<> void Simplex_geom<2>::merge(Nearest_point<2>& nearest, Mat<2, 2> sim, Mat<2> point);
template<> void Simplex_geom<3>::merge(Nearest_point<3>& nearest, Mat<3, 3> sim, Mat<3> point);
//! \cond
template<> void Simplex_geom<3>::visualize(std::string);
//! \endcond

/*! \brief Creates simplices from an ordered array of points representing a polygonal curve.
 * \details `points` must have exactly 2 rows.
 * Each column of `points` is the coordinates of one point.
 * The curve may or may not be closed.
 * Each point represents a unique vertex -- duplicate points are allowed,
 * but you don't need to include each interior point twice.
 * The result can be used to construct a `Simplex_geom<2>`.
 * \see `read_csv()`
 * \relates Simplex_geom
 */
std::vector<Mat<2, 2>> segments(const Mat<dyn, dyn>& points);

}
#endif
