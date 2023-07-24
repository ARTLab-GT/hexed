#ifndef HEXED_NEAREST_POINT_HPP_
#define HEXED_NEAREST_POINT_HPP_

#include "math.hpp"

namespace hexed
{

/*! \brief Helper class for finding the nearest point to a given reference point.
 * \details Implementations of `Surface_geom::nearest_point`
 * often involve searching a set of candidate points
 * to find the one that is nearest to the input point.
 * This class abstracts the logic of comparing the distances and updating the current nearest point
 * to make those implementations more concise.
 */
template <int n_dim = dyn>
class Nearest_point
{
  Mat<n_dim> r;
  Mat<n_dim> p;
  double dist_sq;
  bool e;
  public:
  /*! \brief Initializes an empty point.
   * \details This point is considered "empty" and `point()` won't work
   * until `merge()`ed with a non-empty `Nearest_point`.
   * \param ref Sets `reference()`
   * \param max_distance Optionally initialize the distance to a finite value so that only points within this distance will be considered.
   */
  Nearest_point(Mat<n_dim> ref, double max_distance = huge) : r{ref}, p{ref}, dist_sq{max_distance*max_distance}, e{true} {}
  //! \brief Creates a `Nearest_point` with a specified candidate point and computes the correct distance.
  //! \details This point is non-empty and has a well-defined `point()`.
  Nearest_point(Mat<n_dim> ref, Mat<n_dim> pnt) : r{ref}, p{pnt}, dist_sq{(pnt - ref).squaredNorm()}, e{false} {}
  //! \brief simple copy constructor that just converts the matrices
  template <int m_dim> Nearest_point(const Nearest_point<m_dim>& other)
  {
    if (other.empty()) *this = Nearest_point(other.reference());
    else *this = Nearest_point(other.reference(), other.point());
  }
  Mat<n_dim> reference() const {return r;} //!< reference point that you want to find the nearest point to
  double dist_squared() const {return dist_sq;} //!< squared distance between `reference()` and `point()`
  bool empty() const {return e;} //!< if `true`, `point()` won't work

  //! sets `*this` to the nearest of `*this` and `other`
  void merge(Nearest_point other) {if (other.dist_sq < this->dist_sq) *this = other;}
  //! updates `*this` to point to `new_pnt` if `new_pnt` is closer
  void merge(Mat<n_dim> new_pnt) {merge(Nearest_point(r, new_pnt));}

  //! \brief Current best estimate for nearest point.
  //! \details Throws an exception if `empty()`.
  Mat<n_dim> point() const
  {
    HEXED_ASSERT(!e, "no candidates have been merged");
    return p;
  }
};

}
#endif
