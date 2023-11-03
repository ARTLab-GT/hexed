#ifndef HEXED_SIMPLEX_GEOM_HPP_
#define HEXED_SIMPLEX_GEOM_HPP_

#include "Surface_geom.hpp"
#include "Stopwatch_tree.hpp"
#include "Tree.hpp"

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
  static Mat<n_dim, 2> _get_bounding_box(const std::vector<Mat<n_dim, n_dim>>& sims)
  {
    Mat<n_dim, 2> box;
    box(all, 0).setConstant(huge);
    box(all, 1).setConstant(-huge);
    for (auto& simplex : sims) {
      for (int i_vertex = 0; i_vertex < n_dim; ++i_vertex) {
        box(all, 0) = box(all, 0).cwiseMin(simplex(all, i_vertex));
        box(all, 1) = box(all, 1).cwiseMax(simplex(all, i_vertex));
      }
    }
    return box;
  }

  void sort_simplices(Tree& tree)
  {
    tree.refine();
    std::vector<bool> in_child(tree.misc_data.size(), false);
    for (Tree* child : tree.children()) {
      Mat<n_dim> coords = child->nominal_position();
      for (unsigned i_ind = 0; i_ind < tree.misc_data.size(); ++i_ind) {
        if (!in_child[i_ind]) {
          in_child[i_ind] = true;
          for (int i_vertex = 0; i_vertex < n_dim; ++i_vertex) {
            Mat<n_dim> diff = _simplices[tree.misc_data[i_ind]](all, i_vertex) - coords;
            in_child[i_ind] = in_child[i_ind] && diff.minCoeff() > 0. && diff.maxCoeff() < child->nominal_size();
          }
          if (in_child[i_ind]) child->misc_data.push_back(tree.misc_data[i_ind]);
        }
      }
    }
    std::vector<int> remaining;
    for (unsigned i_ind = 0; i_ind < tree.misc_data.size(); ++i_ind) {
      if (!in_child[i_ind]) remaining.push_back(tree.misc_data[i_ind]);
    }
    tree.misc_data = std::move(remaining);
    for (Tree* child : tree.children()) if (!child->misc_data.empty()) sort_simplices(*child);
  }

  void recursive_nearest(Nearest_point<n_dim>& nearest, Tree& tree, Mat<n_dim> point, double limit)
  {
    if ((point - tree.center()).norm() - tree.nominal_size()*std::sqrt(n_dim)/2 > limit) return;
    for (int i_simplex : tree.misc_data) {
      auto sim = _simplices[i_simplex];
      auto bball = math::bounding_ball(sim);
      if ((bball.center - point).norm() - std::sqrt(bball.radius_sq) <= limit) { // only bother if at least the bounding ball is close enough
        merge(nearest, sim, point);
      }
    }
    for (Tree* child : tree.children()) recursive_nearest(nearest, *child, point, limit);
  }

  void recursive_intersections(std::vector<double>& inters, Tree& tree, Mat<n_dim> point0, Mat<n_dim> point1)
  {
    if (!math::intersects(math::Ball<n_dim>(tree.center(), tree.nominal_size()*n_dim/4.), point0, point1)) return;
    Mat<n_dim> diff = point1 - point0;
    for (int i_simplex : tree.misc_data) {
      auto sim = _simplices[i_simplex];
      if (math::intersects<n_dim>(math::bounding_ball(sim), point0, point1)) { // if the line doesn't even intersect the bounding ball, don't bother
        // compute intersection between line and plane of simplex
        Mat<n_dim, n_dim> lhs;
        lhs(all, 0) = -diff;
        for (int col = 1; col < n_dim; ++col) lhs(all, col) = sim(all, col) - sim(all, 0);
        auto fact = lhs.lu();
        if (fact.rcond() > 1e-7) {
          Mat<n_dim> soln = fact.solve(point0 - sim(all, 0));
          // if intersection is inside simplex, add it to the list
          Eigen::Array<double, n_dim - 1, 1> arr = soln(Eigen::seqN(1, n_dim - 1)).array();
          if ((arr >= -gap_tol).all() && arr.sum() <= 1 + gap_tol) inters.push_back(soln(0));
        }
      }
    }
    for (Tree* child : tree.children()) recursive_intersections(inters, *child, point0, point1);
  }

  std::vector<Mat<n_dim, n_dim>> _simplices;
  Mat<n_dim, 2> _bounding_box;
  Tree _tree;

  public:
  double gap_tol = 1e-6; //!< extend triangles by this amount, relative to their original size, to fill gaps

  Simplex_geom(const std::vector<Mat<n_dim, n_dim>>& sims)
  : _simplices{sims},
    _bounding_box(_get_bounding_box(sims)),
    _tree(n_dim, (_bounding_box(all, 1) - _bounding_box(all, 0)).maxCoeff(), _bounding_box(all, 0))
  {
    for (unsigned i_ind = 0; i_ind < sims.size(); ++i_ind) _tree.misc_data.push_back(i_ind);
    sort_simplices(_tree);
  }

  void visualize(std::string file_name); //!< writes geometry to a Tecplot file with the specified name + file extension

  /*! \details Iterates through all _simplices and finds the nearest point on each,
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
    recursive_nearest(nearest, _tree, point, limit);
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
   * so if the line passes exactly through the shared boundary of multiple _simplices
   * then duplicate intersections may be obtained.
   */
  std::vector<double> intersections(Mat<> point0, Mat<> point1) override
  {
    #if HEXED_OBSESSIVE_TIMING
    Stopwatch sw;
    sw.start();
    #endif
    std::vector<double> inters;
    recursive_intersections(inters, _tree, point0, point1);
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
